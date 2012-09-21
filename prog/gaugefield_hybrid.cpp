#include "gaugefield_hybrid.h"

#include "meta/util.hpp"

void Gaugefield_hybrid::init(int numtasks, cl_device_type primary_device_type)
{

	logger.trace() << "Initialize gaugefield";

	//how many tasks (devices with different purpose) do we need:
	set_num_tasks(numtasks);
	//LZ: for now assume that there is only one device per task

	//allocate memory for private gaugefield on host and initialize (cold start, read in, future: hot start)
	sgf = new Matrixsu3[get_num_gaugefield_elems()];
	init_gaugefield();

	init_devicetypearray(primary_device_type);
	init_opencl();

	this->init_tasks();

	//this has to be done anyways...
	copy_gaugefield_to_all_tasks();
}

void Gaugefield_hybrid::init_devicetypearray(cl_device_type primary_device_type)
{
	//init devicetype array

	//LZ: Note that num_task_types and the input parameter num_dev seem not to fit to each other. However, consider following scenario:
	//    num_task_types is given to this class and can potentially control whether certain additional tasks are performed or not
	//    num_dev controls how many devices should be used/are available (so usually CPU and GPU, but possibly also CPU, GPU from different nodes)
	//    here, the different devices can be assigned automatically to the tasks
	//CP: At the moment, get_num_dev is simply not used, altough it may be in the future...
	//LZ: first application of num_dev: if num_dev==1 then we want to use devices of one type (primary_device_type) only

	devicetypes = new cl_device_type[get_num_tasks()];
	if(get_num_tasks() == 1) {
		devicetypes[0] = primary_device_type;
	} else if (get_num_tasks() == 2) {
		devicetypes[0] = primary_device_type;
		if(primary_device_type == CL_DEVICE_TYPE_GPU) {
			if(get_parameters().get_device_count() == 1)
				devicetypes[1] = CL_DEVICE_TYPE_GPU;
			else
				devicetypes[1] = CL_DEVICE_TYPE_CPU;
		} else {
			if(get_parameters().get_device_count() == 1)
				devicetypes[1] = CL_DEVICE_TYPE_CPU;
			else
				devicetypes[1] = CL_DEVICE_TYPE_GPU;
		}
	} else {
		throw Print_Error_Message("3 or more tasks not yet implemented.");
	}

	logger.debug() << "Wish list for device types:" ;
	for(int n = 0; n < get_num_tasks(); n++) {
		switch(devicetypes[n]) {
			case CL_DEVICE_TYPE_CPU :
				logger.debug() << "CL_DEVICE_TYPE_CPU" ;
				break;
			case CL_DEVICE_TYPE_GPU :
				logger.debug() << "CL_DEVICE_TYPE_GPU" ;
				break;
		}
	}
}

void Gaugefield_hybrid::init_opencl()
{
	cl_int clerr = CL_SUCCESS;

	//Initialize OpenCL,
	logger.trace() << "OpenCL being initialized...";

	cl_uint num_devices_gpu = 0;
	cl_uint num_devices_cpu = 0;
for(auto device: system->get_devices()) {
		switch(device->get_device_type()) {
			case CL_DEVICE_TYPE_GPU:
				++num_devices_gpu;
				break;
			case CL_DEVICE_TYPE_CPU:
				++num_devices_cpu;
				break;
				// ignore other cases
		}
	}
	logger.info() << "\tFound " << num_devices_gpu << " GPU(s) and " << num_devices_cpu << " CPU(s).";

	//LZ: begin debug
	//    num_devices_gpu = 0;
	//  num_devices_cpu = 1;
	//  logger.info() << "Found " << num_devices_gpu << " GPU(s) and " << num_devices_cpu << " CPU(s).";
	// end debug

	set_num_devices(num_devices_gpu + num_devices_cpu);

	//Map tasks and devices
	if(get_num_devices() < 1) throw Print_Error_Message("No devices found.");
	if( num_devices_gpu == 0 ) logger.warn() << "No GPU found.";
	if( num_devices_cpu == 0 ) logger.warn() << "No CPU found.";

	switch (get_num_tasks() ) {

		case 1 :
			if(devicetypes[0] == CL_DEVICE_TYPE_GPU && num_devices_gpu < 1) {
				logger.warn() << "Program demanded GPU device but there is only a CPU device available. Try that." ;
				devicetypes[0] = CL_DEVICE_TYPE_CPU;
			}
			if(devicetypes[0] == CL_DEVICE_TYPE_CPU && num_devices_cpu < 1) {
				logger.warn() << "Program demanded CPU device but there is only a GPU device available. Try that." ;
				devicetypes[0] = CL_DEVICE_TYPE_GPU;
			}
			break;

		case 2 :
			if( get_num_devices() == 1 ) {
				logger.warn() << "You wanted to have two devices, but only one has been found!";
				if( num_devices_gpu == 1 ) {
					devicetypes[0] = CL_DEVICE_TYPE_GPU;
					devicetypes[1] = CL_DEVICE_TYPE_GPU;
				}
				if( num_devices_cpu == 1 ) {
					devicetypes[0] = CL_DEVICE_TYPE_CPU;
					devicetypes[1] = CL_DEVICE_TYPE_CPU;
				}
			}
			break;

		default:
			logger.warn() << "Set of devices could not be mapped properly. Try what happens..." ;
	}

	int len = std::min( get_num_tasks(), get_num_devices() );
	devices                 = new hardware::Device*[len];
	cl_devices              = new cl_device_id[len];

	for(int ntask = 0; ntask < len; ntask++) {
		logger.info() << "\tInitialize device #" << ntask << ":";
		init_devices(ntask);
	}

	//now we need a mapping between devices and tasks for the case of fewer devices than tasks
	device_id_for_task = new int [get_num_tasks()];
	for(int ntask = 0 ; ntask < get_num_tasks() ; ntask++) {
		device_id_for_task[ntask] = ntask;
	}
	if(get_num_devices() < get_num_tasks()) {
		switch (get_num_tasks()) {
			case 2:
				device_id_for_task[1] = 0;
				break;
			default:
				throw Print_Error_Message("Fewer devices than tasks but no proper mapping available.");
		}
	}
}

void Gaugefield_hybrid::init_devices(int ndev)
{
	cl_int clerr = CL_SUCCESS;

	char info[512];

	// currently always use first device of proper type
	// TODO use proper device
	devices[ndev] = 0;
for(auto device: system->get_devices()) {
		if(device->get_device_type() == devicetypes[ndev]) {
			cl_devices[ndev] = device->get_id();
			devices[ndev] = device;
			break;
		}
	}
	if(devices[ndev] == 0) {
		throw Print_Error_Message("Failed to find usable device for task " + ndev);
	}

	clerr = clGetDeviceInfo(cl_devices[ndev], CL_DEVICE_NAME, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\t\tCL_DEVICE_NAME:    " << info;
	clerr = clGetDeviceInfo(cl_devices[ndev], CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.debug() << "\t\t\tCL_DEVICE_VENDOR:  " << info;
	cl_device_type devtype;
	clerr = clGetDeviceInfo(cl_devices[ndev], CL_DEVICE_TYPE, sizeof(cl_device_type), &devtype, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	if(devtype == CL_DEVICE_TYPE_CPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    CPU";
	if(devtype == CL_DEVICE_TYPE_GPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    GPU";
	if(devtype == CL_DEVICE_TYPE_ACCELERATOR) logger.info() << "\t\t\tCL_DEVICE_TYPE:    ACCELERATOR";
	if(devtype != CL_DEVICE_TYPE_CPU && devtype != CL_DEVICE_TYPE_GPU && devtype != CL_DEVICE_TYPE_ACCELERATOR)
		throw Print_Error_Message("Unexpected CL_DEVICE_TYPE...", __FILE__, __LINE__);
	clerr = clGetDeviceInfo(cl_devices[ndev], CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.debug() << "\t\t\tCL_DEVICE_VERSION: " << info;
}



void Gaugefield_hybrid::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
		opencl_modules[ntask] = new Opencl_Module(parameters, get_device_for_task(ntask));
		opencl_modules[ntask]->init();
	}
}


void Gaugefield_hybrid::finalize()
{
	this->finalize_opencl();
	this->delete_variables();
}

void Gaugefield_hybrid::delete_variables()
{
	delete [] device_id_for_task;

	delete [] sgf;

	delete [] devices;
	delete [] cl_devices;

	delete [] devicetypes;

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		delete opencl_modules[ntask];
	}
	delete [] opencl_modules;
}

void Gaugefield_hybrid::finalize_opencl()
{

}

void Gaugefield_hybrid::init_gaugefield()
{
	switch(get_parameters().get_startcondition()) {
		case meta::Inputparameters::start_from_source: {
			sourcefileparameters parameters_source;
			//tmp hmc_gaugefield for filetransfer
			hmc_float * gaugefield_tmp;
			gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float) * NDIM * NC * NC * parameters.get_ntime() * meta::get_volspace(parameters));
			parameters_source.readsourcefile(get_parameters().get_sourcefile().c_str(), get_parameters().get_precision(), &gaugefield_tmp);
			copy_gaugefield_from_ildg_format(get_sgf(), gaugefield_tmp, parameters_source.num_entries_source);
			free(gaugefield_tmp);
		}
		break;
		case meta::Inputparameters::cold_start:
			set_gaugefield_cold(get_sgf());
			break;
		case meta::Inputparameters::hot_start:
			set_gaugefield_hot(get_sgf());
			break;
	}
}

void Gaugefield_hybrid::set_gaugefield_cold(Matrixsu3 * field)
{
	for(int t = 0; t < parameters.get_ntime(); t++) {
		for(size_t n = 0; n < meta::get_volspace(parameters); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const Matrixsu3 tmp = unit_matrixsu3();
				set_to_gaugefield(field, mu, n, t, tmp);
			}
		}
	}
}

void Gaugefield_hybrid::set_gaugefield_hot(Matrixsu3 * field)
{
	for(int t = 0; t < parameters.get_ntime(); t++) {
		for(size_t n = 0; n < meta::get_volspace(parameters); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const Matrixsu3 tmp = random_matrixsu3();
				set_to_gaugefield(field, mu, n, t, tmp);
			}
		}
	}
}

void Gaugefield_hybrid::copy_gaugefield_to_all_tasks()
{
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		copy_gaugefield_to_task(ntask);
	}
}

void Gaugefield_hybrid::copy_gaugefield_to_task(int ntask)
{
	if(ntask < 0 || ntask > get_num_tasks() ) {
		logger.warn() << "Index out of range, copy_gaugefield_to_device does nothing.";
		return;
	}
	opencl_modules[ntask]->importGaugefield(get_sgf());
}

void Gaugefield_hybrid::synchronize(int ntask_reference)
{
	logger.debug() << "Syncrhonizing to data of task " << ntask_reference;
	if(ntask_reference < 0 || ntask_reference > get_num_tasks() ) {
		logger.warn() << "Index out of range, synchronize_gaugefield does nothing.";
		return;
	}
	copy_gaugefield_from_task(ntask_reference);

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		get_device_for_task(ntask)->synchronize();
		if(ntask != ntask_reference) copy_gaugefield_to_task(ntask);
	}
}

void Gaugefield_hybrid::copy_gaugefield_from_task(int ntask)
{
	if(ntask < 0 || ntask > get_num_tasks() ) {
		logger.warn() << "Index out of range, copy_gaugefield_from_device does nothing.";
		return;
	}
	opencl_modules[ntask]->exportGaugefield(get_sgf());
}

hardware::Device* Gaugefield_hybrid::get_device_for_task(int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("task index out of range", __FILE__, __LINE__);
	return devices[device_id_for_task[ntask]];
}

void Gaugefield_hybrid::set_num_tasks (int num)
{
	num_tasks = num;
}

int Gaugefield_hybrid::get_num_tasks ()
{
	return num_tasks;
}

const meta::Inputparameters & Gaugefield_hybrid::get_parameters ()
{
	return  parameters;
}

Matrixsu3 * Gaugefield_hybrid::get_sgf ()
{
	return sgf;
}

void Gaugefield_hybrid::set_sgf (Matrixsu3 * sgf_val)
{
	sgf = sgf_val;
}

void Gaugefield_hybrid::set_num_devices(int num)
{
	num_devices = num;
}

int Gaugefield_hybrid::get_num_devices()
{
	return num_devices;
}

void Gaugefield_hybrid::save(int number)
{
	using namespace std;

	//LZ: generalize the following to larger numbers, if necessary...
	stringstream strnumber;
	strnumber.fill('0');
	strnumber.width(5);
	strnumber << right << number;
	stringstream outfilename;
	outfilename << "conf." << strnumber.str();
	string outputfile = outfilename.str();
	save(outputfile);
}


void Gaugefield_hybrid::save(std::string outputfile)
{
	const size_t NTIME = parameters.get_ntime();
	const size_t gaugefield_buf_size = 2 * NC * NC * NDIM * meta::get_volspace(parameters) * NTIME;
	hmc_float * gaugefield_buf = new hmc_float[gaugefield_buf_size];

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;

	copy_gaugefield_to_ildg_format(gaugefield_buf, get_sgf());

	hmc_float plaq = plaquette();

	int number = 0;

	const size_t NSPACE = parameters.get_nspace();

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NTIME, get_parameters().get_precision(), number, plaq, get_parameters().get_beta(), get_parameters().get_kappa(), get_parameters().get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	delete[] gaugefield_buf;
}


hmc_float Gaugefield_hybrid::plaquette()
{
	hmc_float tdummy;
	hmc_float sdummy;
	//LZ: calculation of tdummy and sdummy is unnecessary for this function, chance to speed up a little bit...
	return plaquette(&tdummy, &sdummy);
}

void Gaugefield_hybrid::print_gaugeobservables(int iter)
{
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	hmc_complex pol = polyakov();
	logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
}

void Gaugefield_hybrid::print_gaugeobservables(int iter, std::string filename)
{
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	hmc_complex pol = polyakov();
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) throw File_Exception(filename);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	hmc_float plaq  = 0;
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_complex pol;
	opencl_modules[ntask]->gaugeobservables(&plaq, &tplaq, &splaq, &pol);
	logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask, std::string filename)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	hmc_float plaq  = 0;
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_complex pol;
	opencl_modules[ntask]->gaugeobservables(&plaq, &tplaq, &splaq, &pol);
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) throw File_Exception(filename);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
}

#ifdef _PROFILING_
void Gaugefield_hybrid::print_profiling(std::string filename)
{
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
		opencl_modules[ntask]->print_profiling(filename, ntask);
	}
}
#endif

void Gaugefield_hybrid::set_to_gaugefield(Matrixsu3 * field, const size_t mu, const size_t x, const size_t t, const Matrixsu3 val)
{
	size_t pos = get_global_link_pos(mu, x, t, parameters);
	field[pos] = val;
}

Matrixsu3 Gaugefield_hybrid::get_from_gaugefield(const Matrixsu3 * field, const size_t mu, const size_t x, const size_t t) const
{
	return field[get_global_link_pos(mu, x, t, parameters)];
}

size_t Gaugefield_hybrid::get_num_gaugefield_elems() const
{
	return NDIM * meta::get_volspace(parameters) * parameters.get_ntime();
}

void Gaugefield_hybrid::copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, hmc_float * gaugefield_tmp, int check)
{
	//little check if arrays are big enough
	if (meta::get_vol4d(parameters) *NDIM * NC * NC * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t NSPACE = parameters.get_nspace();
	int cter = 0;
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in a complex array
						hmc_complex tmp [NC][NC];
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								int pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								tmp[m][n].re = gaugefield_tmp[pos];
								tmp[m][n].im = gaugefield_tmp[pos + 1];
								cter++;
							}
						}
						//store su3matrix tmp in our format
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);

						//copy hmc_su3matrix to Matrixsu3 format
						Matrixsu3 destElem;
						destElem.e00 = tmp[0][0];
						destElem.e01 = tmp[0][1];
						destElem.e02 = tmp[0][2];
						destElem.e10 = tmp[1][0];
						destElem.e11 = tmp[1][1];
						destElem.e12 = tmp[1][2];
						destElem.e20 = tmp[2][0];
						destElem.e21 = tmp[2][1];
						destElem.e22 = tmp[2][2];

						set_to_gaugefield(gaugefield, (l + 1) % NDIM, spacepos, t, destElem);
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	return;
}

void Gaugefield_hybrid::copy_gaugefield_to_ildg_format(hmc_float * dest, Matrixsu3 * source_in)
{
	const size_t NSPACE = parameters.get_nspace();
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);
						hmc_complex destElem [NC][NC];

						Matrixsu3 srcElem = get_from_gaugefield(source_in, (l + 1) % NDIM, spacepos, t);
						destElem[0][0] = srcElem.e00;
						destElem[0][1] = srcElem.e01;
						destElem[0][2] = srcElem.e02;
						destElem[1][0] = srcElem.e10;
						destElem[1][1] = srcElem.e11;
						destElem[1][2] = srcElem.e12;
						destElem[2][0] = srcElem.e20;
						destElem[2][1] = srcElem.e21;
						destElem[2][2] = srcElem.e22;

						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								dest[pos]     = destElem[m][n].re;
								dest[pos + 1] = destElem[m][n].im;
							}
						}
					}
				}
			}
		}
	}

	return;
}

//taken from the opencl files
hmc_complex trace_matrixsu3(const Matrixsu3 p)
{
	hmc_complex out;
	out.re = p.e00.re;
	out.im = p.e00.im;
	out.re += p.e11.re;
	out.im += p.e11.im;
	out.re += p.e22.re;
	out.im += p.e22.im;
	return out;
}

hmc_float Gaugefield_hybrid::plaquette(hmc_float* tplaq, hmc_float* splaq)
{
	hmc_float plaq = 0;
	*tplaq = 0;
	*splaq = 0;
	int coord[NDIM];

	for(int t = 0; t < parameters.get_ntime(); t++) {
		for(int x = 0; x < parameters.get_nspace(); x++) {
			for(int y = 0; y < parameters.get_nspace(); y++) {
				for(int z = 0; z < parameters.get_nspace(); z++) {
					for(int mu = 0; mu < NDIM; mu++) {
						for(int nu = 0; nu < mu; nu++) {
							coord[0] = t;
							coord[1] = x;
							coord[2] = y;
							coord[3] = z;
							Matrixsu3 prod;
							prod = local_plaquette(get_sgf(), coord, mu, nu, parameters );
							hmc_float tmpfloat = trace_matrixsu3(prod).re;
							plaq += tmpfloat;
							if(mu == 0 || nu == 0) {
								*tplaq += tmpfloat;
							} else {
								*splaq += tmpfloat;
							}
						}
					}
				}
			}
		}
	}

	*tplaq /= static_cast<hmc_float>(meta::get_vol4d(parameters) * NC * (NDIM - 1));
	*splaq /= static_cast<hmc_float>(meta::get_vol4d(parameters) * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	return plaq * 2.0 / static_cast<hmc_float>(meta::get_vol4d(parameters) * NDIM * (NDIM - 1) * NC);
}

hmc_complex Gaugefield_hybrid::polyakov()
{
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int n = 0; n < meta::get_volspace(parameters); n++) {
		Matrixsu3 prod;
		prod = local_polyakov(get_sgf(), n, parameters);
		hmc_complex tmpcomplex = trace_matrixsu3(prod);
		res.re += tmpcomplex.re;
		res.im += tmpcomplex.im;
	}

	res.re /= static_cast<hmc_float>(NC * meta::get_volspace(parameters));
	res.im /= static_cast<hmc_float>(NC * meta::get_volspace(parameters));
	return res;
}

