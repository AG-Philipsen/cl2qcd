#include "opencl_fermions.h"
#include "logger.hpp" 

hmc_error Opencl_fermions::fill_kernels_file ()
{
	//give a list of all kernel-files
	Opencl::fill_kernels_file();

	cl_kernels_file.push_back("types_fermions.h");
	cl_kernels_file.push_back("operations_su3vec.cl");
	cl_kernels_file.push_back("operations_spinor.cl");
	cl_kernels_file.push_back("operations_spinorfield.cl");
	#ifdef _USE_GPU_
	cl_kernels_file.push_back("operations_fermionmatrix_GPU.cl");
	#else
	cl_kernels_file.push_back("operations_fermionmatrix.cl");
	#endif
	if(get_parameters()->get_use_eo() == TRUE){
		cl_kernels_file.push_back("operations_spinorfield_eo.cl");
		#ifdef _USE_GPU_
		cl_kernels_file.push_back("operations_fermionmatrix_eo_GPU.cl");
		#else
		cl_kernels_file.push_back("operations_fermionmatrix_eo.cl");
		#endif
	}
	cl_kernels_file.push_back("fermionobservables.cl");	

	return HMC_SUCCESS;  
}

hmc_error Opencl_fermions::fill_collect_options(stringstream* collect_options)
{

	Opencl::fill_collect_options(collect_options);
	*collect_options << " -D_FERMIONS_" << " -DSPINORSIZE=" << SPINORSIZE << " -DHALFSPINORSIZE=" << HALFSPINORSIZE << " -DSPINORFIELDSIZE=" << SPINORFIELDSIZE << " -DEOPREC_SPINORFIELDSIZE=" << EOPREC_SPINORFIELDSIZE;


	switch (get_parameters()->get_fermact()) {
		case TWISTEDMASS :
			*collect_options << " -D_TWISTEDMASS_";
			break;
		case CLOVER :
			*collect_options << " -D_CLOVER_";
			break;
	}

	//CP: give kappa and its negative value
	hmc_float kappa_tmp = get_parameters()->get_kappa();
	*collect_options << " -DKAPPA=" << kappa_tmp;
	*collect_options << " -DMKAPPA=" << -kappa_tmp;
	//CP: These are the kappas including BC in spatial and temporal direction
	hmc_float tmp_spatial = (get_parameters()->get_theta_fermion_spatial() * PI)/( (hmc_float) NSPACE);
	hmc_float tmp_temporal = (get_parameters()->get_theta_fermion_temporal() * PI)/( (hmc_float) NTIME);
	//BC: on the corners in each direction: exp(i theta) -> on each site exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)
	*collect_options << " -DKAPPA_SPATIAL_RE=" << kappa_tmp*cos(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_RE=" << -kappa_tmp*cos(tmp_spatial);	
	*collect_options << " -DKAPPA_SPATIAL_IM=" << kappa_tmp*sin(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_IM=" << -kappa_tmp*sin(tmp_spatial);		
	
	*collect_options << " -DKAPPA_TEMPORAL_RE=" << kappa_tmp*cos(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_RE=" << -kappa_tmp*cos(tmp_temporal);	
	*collect_options << " -DKAPPA_TEMPORAL_IM=" << kappa_tmp*sin(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_IM=" << -kappa_tmp*sin(tmp_temporal);		
	
	switch (get_parameters()->get_fermact()) {
	case TWISTEDMASS :
		*collect_options << " -DMU=" << get_parameters()->get_mu();
	  get_parameters()->calc_mubar();
	  *collect_options << " -DMUBAR=" << get_parameters()->get_mubar();
	  get_parameters()->set_mubar_negative();
	  *collect_options << " -DMMUBAR=" << get_parameters()->get_mubar();
	  get_parameters()->set_mubar_negative();
	  break;
	case CLOVER :
	  *collect_options << " -DCSW=" << get_parameters()->get_csw();
	  break;
	}

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::fill_buffers()
{
	Opencl::fill_buffers();

	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sense on CPU

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups;

	logger.trace() << "init buffer for solver...";
	int clerr = CL_SUCCESS;

	//	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	//int eoprec_spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
	int eoprec_spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
	int global_buf_size = complex_size * num_groups;
	int global_buf_size_float = float_size * num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;

  clmem_corr = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_corr failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_inout = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_source = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_source failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rn = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rn failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rhat = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rhat failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_v = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_p = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_p failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_s = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_s failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_t = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_t failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_aux = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //LZ only use the following if we want to apply even odd preconditioning
  if(get_parameters()->get_use_eo()==TRUE) {
    clmem_inout_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_source_even = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_source_odd = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }  
    clmem_rn_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_rhat_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_v_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_p_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_s_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_t_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }  
    clmem_aux_eoprec = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_tmp_eoprec_1 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }  
    clmem_tmp_eoprec_2 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clmem_tmp_eoprec_3 = clCreateBuffer(context,CL_MEM_READ_WRITE,eoprec_spinorfield_size,0,&clerr);;
    if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_inout_eoprec failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    } 

  } //end if: eoprec
  
  clmem_rho = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
      cout<<"creating clmem_rho failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
  }
  clmem_rho_next = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rho_next failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_alpha = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_alpha failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_omega = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_omega failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_beta = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_beta failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_tmp1 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_tmp2 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_one = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_one failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_minusone = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_minusone failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

  clmem_scalar_product_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_scalar_product_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
    
  clmem_resid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_resid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
	clmem_trueresid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_trueresid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_global_squarenorm_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size_float,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_global_squarenorm_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  

  clerr = clEnqueueWriteBuffer(queue,clmem_one,CL_TRUE,0,complex_size,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_one failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  } 
  clerr = clEnqueueWriteBuffer(queue,clmem_minusone,CL_TRUE,0,complex_size,&minusone,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_minusone failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::fill_kernels()
{
	int clerr = HMC_SUCCESS;
	
	Opencl::fill_kernels();
	
	logger.debug() << "Create fermion kernels...";
	Qplus = clCreateKernel(clprogram,"Qplus",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating Qplus kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( Qplus );
	Qminus = clCreateKernel(clprogram,"Qminus",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating Qminus kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( Qminus );
	gamma5 = clCreateKernel(clprogram,"gamma5",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating gamma5 kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	if( logger.beDebug() )
		printResourceRequirements( gamma5 );
	}
	M = clCreateKernel(clprogram,"M",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating M kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( M );
	ps_correlator = clCreateKernel(clprogram,"ps_correlator",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating ps_correlator kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( ps_correlator );
	set_spinorfield_cold = clCreateKernel(clprogram,"set_spinorfield_cold",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating set_spinorfield_cold kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( set_spinorfield_cold );
	saxpy = clCreateKernel(clprogram, "saxpy", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating saxpy kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( saxpy );
	saxsbypz = clCreateKernel(clprogram, "saxsbypz", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating saxsbypz kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( saxsbypz );
	scalar_product = clCreateKernel(clprogram, "scalar_product", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating scalar_product kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( scalar_product );
	scalar_product_reduction = clCreateKernel(clprogram, "scalar_product_reduction", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating scalar_product_reduction kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( scalar_product_reduction );
	set_zero_spinorfield = clCreateKernel(clprogram, "set_zero_spinorfield", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating set_zero_spinorfield kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( set_zero_spinorfield );
	global_squarenorm = clCreateKernel(clprogram, "global_squarenorm", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating global_squarenorm kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( global_squarenorm );
	global_squarenorm_reduction = clCreateKernel(clprogram, "global_squarenorm_reduction", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating global_squarenorm_reduction kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( global_squarenorm_reduction );
	ratio = clCreateKernel(clprogram, "ratio", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating ratio kernel failed, aborting. " << clerr << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( ratio );
	product = clCreateKernel(clprogram, "product", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating product kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( product );
	convert_to_kappa_format = clCreateKernel(clprogram, "convert_to_kappa_format", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating convert_to_kappa_format kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( convert_to_kappa_format );
	convert_from_kappa_format = clCreateKernel(clprogram, "convert_from_kappa_format", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating convert_from_kappa_format kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( convert_from_kappa_format );
	create_point_source = clCreateKernel(clprogram, "create_point_source", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating create_point_source kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( create_point_source );

	//Kernels needed if eoprec is used
	if(get_parameters()->get_use_eo() == TRUE) {
		M_sitediagonal = clCreateKernel(clprogram, "M_sitediagonal", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating M_sitediagonal kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( M_sitediagonal );
		M_inverse_sitediagonal = clCreateKernel(clprogram, "M_inverse_sitediagonal", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating M_inverse_sitediagonal kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
			printResourceRequirements( M_inverse_sitediagonal );		
		convert_from_eoprec = clCreateKernel(clprogram,"convert_from_eoprec",&clerr);
		if(clerr!=CL_SUCCESS) {
			cout<<"...creating convert_from_eoprec kernel failed, aborting."<<endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( convert_from_eoprec );
		set_eoprec_spinorfield_cold = clCreateKernel(clprogram,"set_eoprec_spinorfield_cold",&clerr);
		if(clerr!=CL_SUCCESS) {
			cout<<"...creating set_eoprec_spinorfield_cold kernel failed, aborting."<<endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( set_eoprec_spinorfield_cold );
		gamma5_eoprec = clCreateKernel(clprogram,"gamma5_eoprec",&clerr);	
		if(clerr!=CL_SUCCESS) {
			cout<<"...creating gamma5_eoprec kernel failed, aborting."<<endl;
			exit(HMC_OCLERROR);
		}	
		if( logger.beDebug() )
			printResourceRequirements( gamma5_eoprec );
		convert_to_kappa_format_eoprec = clCreateKernel(clprogram, "convert_to_kappa_format_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating convert_to_kappa_format_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( convert_to_kappa_format_eoprec );
		convert_from_kappa_format_eoprec = clCreateKernel(clprogram, "convert_from_kappa_format_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating convert_from_kappa_format_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( convert_from_kappa_format_eoprec );
		dslash_eoprec = clCreateKernel(clprogram, "dslash_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating dslash_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( dslash_eoprec );
		saxpy_eoprec = clCreateKernel(clprogram, "saxpy_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating saxpy_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( saxpy_eoprec );
		saxsbypz_eoprec = clCreateKernel(clprogram, "saxsbypz_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating saxsbypz_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( saxsbypz_eoprec );
		scalar_product_eoprec = clCreateKernel(clprogram, "scalar_product_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating scalar_product_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( scalar_product_eoprec );
		set_zero_spinorfield_eoprec = clCreateKernel(clprogram, "set_zero_spinorfield_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating set_zero_spinorfield_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( set_zero_spinorfield_eoprec );
		global_squarenorm_eoprec = clCreateKernel(clprogram, "global_squarenorm_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating global_squarenorm_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
		if( logger.beDebug() )
		printResourceRequirements( global_squarenorm_eoprec );
		create_point_source_eoprec = clCreateKernel(clprogram, "create_point_source_eoprec", &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "...creating create_point_source_eoprec kernel failed, aborting." << endl;
			exit(HMC_OCLERROR);
		}
	}
	
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters)
{
	(*timer).reset();
	hmc_error err = Opencl::init(wanted_device_type, timer, parameters);
	return err;
	(*timer).add();
}


hmc_error Opencl_fermions::convert_to_kappa_format_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(convert_to_kappa_format, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	clerr = clEnqueueNDRangeKernel(queue, convert_to_kappa_format, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue convert_to_kappa_format kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(convert_from_kappa_format, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(convert_from_kappa_format, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	clerr = clEnqueueNDRangeKernel(queue, convert_from_kappa_format, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue convert_from_kappa_format kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_from_eoprec,0,sizeof(cl_mem),&in1); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(convert_from_eoprec,1,sizeof(cl_mem),&in2); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(convert_from_eoprec,2,sizeof(cl_mem),&out); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 
  clerr = clEnqueueNDRangeKernel(queue, convert_from_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_from_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 clFinish(queue);
  
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error Opencl_fermions::convert_to_kappa_format_eoprec_device(cl_mem in, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(convert_to_kappa_format_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 
  clerr = clEnqueueNDRangeKernel(queue, convert_to_kappa_format_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue convert_to_kappa_format_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(convert_from_kappa_format_eoprec, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(convert_from_kappa_format_eoprec, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	clerr = clEnqueueNDRangeKernel(queue, convert_from_kappa_format_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue convert_from_kappa_format_eoprec kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_spinorfield_to_device(spinorfield* host_spinorfield,  usetimer* timer){
  (*timer).reset();

  /** @todo: spinorfield_size should propably be private */
	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	//int clerr = clEnqueueWriteBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
	int clerr = clEnqueueWriteBuffer(queue,clmem_corr,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_eoprec_spinorfield_to_device(spinorfield_eoprec* host_spinorfield,  usetimer* timer){
  (*timer).reset();

	int spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_inout_eoprec,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_source_to_device(spinorfield* host_source,  usetimer* timer){
  (*timer).reset();

	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_source,CL_TRUE,0,spinorfield_size,host_source,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_eoprec_source_to_device(spinorfield_eoprec* host_source1, spinorfield_eoprec* host_source2, usetimer* timer){
  (*timer).reset();

	int spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_source_even,CL_TRUE,0,spinorfield_size,host_source1,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_source_odd,CL_TRUE,0,spinorfield_size,host_source2,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
		
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::get_spinorfield_from_device(spinorfield* host_spinorfield, usetimer* timer){
  (*timer).reset();

	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
//   int clerr = clEnqueueReadBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
	int clerr = clEnqueueReadBuffer(queue,clmem_corr,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error Opencl_fermions::get_eoprec_spinorfield_from_device(spinorfield_eoprec* host_spinorfield, usetimer* timer){
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*EOPREC_SPINORFIELDSIZE;
  int clerr = clEnqueueReadBuffer(queue,clmem_inout_eoprec,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, spinorfield_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_eoprec_spinor_device(cl_mem in, cl_mem out, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, spinorfield_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	hmc_float tmp;
	clerr = clEnqueueReadBuffer(queue, in, CL_TRUE, 0, sizeof(hmc_float), &tmp, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		cout << "errorcode :" << clerr << endl;
		exit(HMC_OCLERROR);
	}
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	hmc_complex tmp;
	clerr = clEnqueueReadBuffer(queue, in, CL_TRUE, 0, sizeof(hmc_complex), &tmp, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		cout << "errorcode :" << clerr << endl;
		exit(HMC_OCLERROR);
	}
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::Qplus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  (*timer).reset();
  int clerr =CL_SUCCESS;

  clerr = clSetKernelArg(Qplus,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(Qplus,1,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(Qplus,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,Qplus,1,0,&global_work_size,&local_work_size,0,0,NULL);

  (*timer).add();

   clFinish(queue);

  return HMC_SUCCESS;

}

hmc_error Opencl_fermions::Qminus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  (*timer).reset();
  int clerr =CL_SUCCESS;

  clerr = clSetKernelArg(Qminus,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(Qminus,1,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(Qminus,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,Qminus,1,0,&global_work_size,&local_work_size,0,0,NULL);

  (*timer).add();

   clFinish(queue);

  return HMC_SUCCESS;

}

hmc_error Opencl_fermions::QplusQminus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  int clerr =CL_SUCCESS;

	/** @todo one could save one field here if an additional copying would be included in the end... */
  clerr = Qminus_device(in, clmem_tmp, local_work_size, global_work_size, timer);
  if(clerr!=CL_SUCCESS) {
    cout<<"Qminus_device failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = Qplus_device(clmem_tmp, out, local_work_size, global_work_size, timer);
  if(clerr!=CL_SUCCESS) {
    cout<<"Qplus_device failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
  return HMC_SUCCESS;

}

hmc_error Opencl_fermions::M_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer* dslashtimer, usetimer* Mdiagtimer){
// 	(*timer).reset();
	
  int clerr =CL_SUCCESS;
  cl_event event;
  clerr = clSetKernelArg(M,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M,1,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M,1,0,&global_work_size,&local_work_size,0,0,&event);

	int done = clWaitForEvents(1, &event);

	//CP: Old method
//   clFinish(queue);
//   (*timer).add();
	
// 	cout << "wait for events gave: " << done << " and it should be " << CL_COMPLETE << endl;
	//Get Event-Infos
  cl_ulong time_start;
  cl_ulong time_end;
  cl_ulong time_init;
  cl_ulong time_queue;
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &time_start, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &time_end, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &time_queue, NULL);
  clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &time_init, NULL);

  cout.precision(64);
  
//   cout << "times measured by clProfiling"<< endl;
//   cout << "Queue Time: " << time_queue << endl;
//   cout << "Submit Time: " << time_init << endl;
//   cout << "Start Time: " << time_start << endl;
//   cout << "End Time: " << time_end  << endl;
//   cout << "Resulting total Time: " << (time_end - time_queue )*0.001<< endl;
// 	cout << "Resulting Execution Time: " << (time_end - time_start )*0.001<< endl;
// 	cout << "time measured by use_timer:  " << noop.getTime() << endl;  

	//calc Execution Time, in microsecs
	uint64_t exec_time = (time_end - time_start )*0.001;
	(*timer).add(exec_time);

  return HMC_SUCCESS;

}


hmc_error Opencl_fermions::gamma5_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer *timer){
	(*timer).reset();
	int clerr =CL_SUCCESS;

  clerr = clSetKernelArg(gamma5,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(gamma5,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,gamma5,1,0,&global_work_size,&local_work_size,0,0,NULL);
	
	(*timer).add();
}

hmc_error Opencl_fermions::gamma5_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer *timer){
	(*timer).reset();
	int clerr =CL_SUCCESS;

  clerr = clSetKernelArg(gamma5_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(gamma5_eoprec,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,gamma5_eoprec,1,0,&global_work_size,&local_work_size,0,0,NULL);
	
	(*timer).add();
}

hmc_error Opencl_fermions::Aee_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer * singletimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * latimer)
{
	(*timer).reset();

	int even = EVEN;
	int odd = ODD;
	const size_t ls = local_work_size;
	const size_t gs = global_work_size;

	dslash_eoprec_device(in, clmem_tmp_eoprec_1, odd, ls, gs, dslashtimer);
	M_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, ls, gs, Mdiagtimer);
	dslash_eoprec_device(clmem_tmp_eoprec_2, out, even, ls, gs, dslashtimer);
	M_sitediagonal_device(in, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);

	copy_eoprec_spinor_device(out, clmem_tmp_eoprec_3, singletimer);

	saxpy_eoprec_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1, clmem_one, out, ls, gs, latimer);

//  int clerr =CL_SUCCESS;
//  clerr = clSetKernelArg(saxpy_eoprec,0,sizeof(cl_mem),&clmem_tmp_eoprec_3);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,1,sizeof(cl_mem),&clmem_tmp_eoprec_1);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,2,sizeof(cl_mem),&clmem_one);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clSetKernelArg(saxpy_eoprec,3,sizeof(cl_mem),&out);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clerr = clEnqueueNDRangeKernel(queue,saxpy_eoprec,1,0,&gs,&ls,0,0,NULL);
//   if(clerr!=CL_SUCCESS) {
//     cout<<"enqueue saxpy_eoprec kernel failed, aborting..."<<endl;
//     exit(HMC_OCLERROR);
//   }
//   clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::dslash_eoprec_device(cl_mem in, cl_mem out, int evenodd, const size_t local_work_size, const size_t global_work_size, usetimer * timer)
{
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;

	int eo = evenodd;

	clerr = clSetKernelArg(dslash_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,2,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,3,sizeof(int),&eo);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,dslash_eoprec,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
   clFinish(queue);	
	(*timer).add();

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::M_inverse_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer)
{

	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(M_inverse_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(M_inverse_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, M_inverse_sitediagonal, 1, 0, &gs, &ls, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue M_inverse_sitediagonal kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::M_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer)
{

	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(M_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(M_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, M_sitediagonal, 1, 0, &gs, &ls, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue M_sitediagonal kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(saxpy, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, saxpy, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue saxpy kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_spinorfield_cold_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(set_spinorfield_cold,0,sizeof(cl_mem),&clmem_inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,set_spinorfield_cold,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue ps_correlator kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
   clFinish(queue);
	
	(*timer).add();
	
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_eoprec_spinorfield_cold_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(set_eoprec_spinorfield_cold,0,sizeof(cl_mem),&clmem_inout_eoprec); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,set_eoprec_spinorfield_cold,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue ps_correlator kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 clFinish(queue);
	
	(*timer).add();
	
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(saxpy_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, saxpy_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue saxpy_eoprec kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(saxsbypz, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, saxsbypz, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue saxsbypz kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(saxsbypz_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, saxsbypz_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue saxsbypz_eoprec kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(scalar_product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(scalar_product, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product, 3, sizeof(hmc_complex) * local_work_size, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, scalar_product, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);
	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, scalar_product_reduction, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product_reduction kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(scalar_product_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 3, sizeof(hmc_complex) * local_work_size, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, scalar_product_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);
	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, scalar_product_reduction, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product_reduction kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(ratio, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(ratio, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(ratio, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP:this needs only one kernel!!
	size_t one = 1;
	clerr = clEnqueueNDRangeKernel(queue, ratio, 1, 0, &one, &one, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue ratio kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(product, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP:this needs only one kernel!!
	size_t one = 1;
	clerr = clEnqueueNDRangeKernel(queue, product, 1, 0, &one, &one, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue product kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(global_squarenorm, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm, 2, sizeof(hmc_float) * local_work_size, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, global_squarenorm, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue global_squarenorm kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);
	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, global_squarenorm_reduction, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product_reduction kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();

	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(global_squarenorm_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_eoprec, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_eoprec, 2, sizeof(hmc_float) * local_work_size, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, global_squarenorm_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue global_squarenorm kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);
	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, global_squarenorm_reduction, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue scalar_product_reduction kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_zero_spinorfield_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(set_zero_spinorfield, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, set_zero_spinorfield, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue set_zero_spinorfield kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::set_zero_spinorfield_eoprec_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clSetKernelArg(set_zero_spinorfield_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, set_zero_spinorfield_eoprec, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue set_zero_spinorfield_eoprec kernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::copy_complex_device(cl_mem in, cl_mem out, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int complex_size = sizeof(hmc_complex);

	clerr = clEnqueueCopyBuffer(queue, in, out, 0, 0, complex_size , 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax)
{

	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter=0; iter<cgmax; iter++){
	        if(iter%iter_refresh==0) {
			set_zero_spinorfield_device(clmem_v, localsize, globalsize, latimer); 
			set_zero_spinorfield_device(clmem_p, localsize, globalsize, latimer);
			M_device(clmem_inout, clmem_rn, localsize, globalsize, Mtimer, dslashtimer, Mdiagtimer);

			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
			copy_spinor_device(clmem_rn, clmem_rhat, singletimer);

			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			copy_complex_device(clmem_one, clmem_omega, singletimer);
			copy_complex_device(clmem_one, clmem_rho, singletimer);

			//CP: calc initial residuum for output, this is not needed for the algorithm!!
 			//set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
 			//copy_float_from_device(clmem_resid, &resid, copytimer);
			//cout << "initial residuum at iter " << iter << "is: " << scientific << resid << endl;
			//printf("initial residuum at iter %i is %.40e\n", iter, resid);
		}

		set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
		copy_complex_device(clmem_rho_next, clmem_rho, singletimer);
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1, singletimer);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2, singletimer);
		saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p, local_work_size, global_work_size, latimer);

		M_device(clmem_p, clmem_v, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);

		set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha, singletimer);

		saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s, local_work_size, global_work_size, latimer);

		M_device(clmem_s, clmem_t, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);

		set_complex_to_scalar_product_device(clmem_t, clmem_s, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		//!!CP: this can also be global_squarenorm
		set_complex_to_scalar_product_device(clmem_t, clmem_t, clmem_tmp2, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega, singletimer);

		saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn, local_work_size, global_work_size, latimer);

		saxsbypz_device(clmem_p, clmem_s, clmem_inout, clmem_alpha, clmem_omega, clmem_inout, local_work_size, global_work_size, latimer);

		set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid<epssquare) {	
		  M_device(clmem_inout,clmem_aux,local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);
			saxpy_device(clmem_aux, clmem_source, clmem_one, clmem_aux, local_work_size, global_work_size, latimer); 
			set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_trueresid, &trueresid, copytimer);
 			//cout << "\tsolver converged! residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid<epssquare)
				return HMC_SUCCESS;
		}
		else{
		  //printf("residuum at iter%i is:\t%.10e\n", iter, resid);//cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::bicgstab_eoprec_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax)
{
	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter = 0; iter < cgmax; iter++) {

		if(iter % iter_refresh == 0) {
			set_zero_spinorfield_eoprec_device(clmem_v_eoprec, localsize, globalsize, latimer);
			set_zero_spinorfield_eoprec_device(clmem_p_eoprec, localsize, globalsize, latimer);

			Aee_device(clmem_inout_eoprec, clmem_rn_eoprec, localsize, globalsize, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);

			saxpy_eoprec_device(clmem_rn_eoprec, clmem_source_even, clmem_one, clmem_rn_eoprec, localsize, globalsize, latimer);
			copy_eoprec_spinor_device(clmem_rn_eoprec, clmem_rhat_eoprec, singletimer);

			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			copy_complex_device(clmem_one, clmem_omega, singletimer);
			copy_complex_device(clmem_one, clmem_rho, singletimer);

			//!!CP: calc initial residuum, this is not needed for the algorithm!!
 			//set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
 			//copy_float_from_device(clmem_resid, &resid, copytimer);
 			//cout << "initial residuum is: " << resid << endl;
		}

		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
		copy_complex_device(clmem_rho_next, clmem_rho, singletimer);
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1, singletimer);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2, singletimer);
		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec, local_work_size, global_work_size, latimer);

		Aee_device(clmem_p_eoprec, clmem_v_eoprec, local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);

		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_v_eoprec, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha, singletimer);

		saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_s_eoprec, local_work_size, global_work_size, latimer);

		Aee_device(clmem_s_eoprec, clmem_t_eoprec, local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);

		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		//!!CP: can this also be global_squarenorm??
		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_t_eoprec, clmem_tmp2, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega, singletimer);

		saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec, local_work_size, global_work_size, latimer);

		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_s_eoprec, clmem_inout_eoprec, clmem_alpha, clmem_omega, clmem_inout_eoprec, local_work_size, global_work_size, latimer);

		set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid < epssquare) {
			Aee_device(clmem_inout_eoprec, clmem_aux_eoprec, local_work_size, global_work_size, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			saxpy_eoprec_device(clmem_aux_eoprec, clmem_source_even, clmem_one, clmem_aux_eoprec, local_work_size, global_work_size, latimer);
			set_float_to_global_squarenorm_eoprec_device(clmem_aux_eoprec, clmem_trueresid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_trueresid, &trueresid, copytimer);
 			//cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid<epssquare)
				return HMC_SUCCESS;
		}
		else{
		  // 			cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax)
{
	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	int iter;
	for(iter = 0; iter < cgmax; iter ++) {
		if(iter % iter_refresh == 0) {
			M_device(clmem_inout, clmem_rn, localsize, globalsize, Mtimer, dslashtimer, Mdiagtimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
			copy_spinor_device(clmem_rn, clmem_p, copytimer);
		}
		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_device(clmem_p, clmem_p, clmem_omega, local_work_size, global_work_size, scalarprodtimer);
		//A pn --> v
		M_device(clmem_p, clmem_v, local_work_size, global_work_size, Mtimer, dslashtimer, Mdiagtimer);
		set_complex_to_scalar_product_device(clmem_p, clmem_v, clmem_rho, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha, singletimer);

		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1, singletimer);

		saxpy_device(clmem_inout, clmem_p, clmem_tmp1, clmem_inout, localsize, globalsize, latimer);
		//rn+1 -> rhat
		saxpy_device(clmem_rn, clmem_v, clmem_alpha, clmem_rhat, localsize, globalsize, latimer);

		set_float_to_global_squarenorm_device(clmem_rhat, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid < epssquare) {
			return HMC_SUCCESS;
		} else {
			//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
			set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta, singletimer);

			//pn+1 = rn+1 + beta*pn
			set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2, singletimer);
			saxpy_device(clmem_p, clmem_rhat, clmem_tmp2, clmem_p, localsize, globalsize, latimer);
		}
	}
	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::solver_device(usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax){
	(*solvertimer).reset();
	convert_to_kappa_format_device(clmem_inout, ls, gs, latimer);
	bicgstab_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer , ls, gs, cgmax);
	convert_from_kappa_format_device(clmem_inout, clmem_inout, ls, gs, latimer);
	(*solvertimer).add();

	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::solver_eoprec_device(usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax){
	(*solvertimer).reset();
	
	//CP: even solution
	convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, ls, gs, latimer);
	bicgstab_eoprec_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer ,dslashtimer, Mdiagtimer, ls, gs, cgmax);	

	//P: odd solution

	/** @todo CP: perhaps one can save some variables used here */
	dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_3, ODD, ls, gs, dslashtimer);
	M_inverse_sitediagonal_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);
	M_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_2, ls, gs, Mdiagtimer);
	saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, clmem_one, clmem_tmp_eoprec_3,  ls, gs, latimer);

	convert_from_kappa_format_eoprec_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_3,  ls, gs, latimer);
	convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, ls, gs, latimer);
	//CP: whole solution
	//CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eoprec_3
	usetimer noop;
	convert_from_eoprec_device(clmem_inout_eoprec,clmem_tmp_eoprec_3,clmem_inout, ls, gs, &noop);
	
	(*solvertimer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::create_point_source_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer)
{

	set_zero_spinorfield_device(clmem_source, ls, gs, latimer);

	(*latimer).reset();
	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), &clmem_source);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueNDRangeKernel(queue, create_point_source, 1, 0, &gs, &ls, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "enqueue create_point_source kernel failed (" << clerr << "), aborting..." << endl;
		exit(HMC_OCLERROR);
	}
 	clFinish(queue);

	(*latimer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::create_point_source_eoprec_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer)
{

	set_zero_spinorfield_eoprec_device(clmem_source_even, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_source_odd, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_2, ls, gs, latimer);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_1, ls, gs, latimer);

	int glob_pos = get_global_pos(spacepos, timepos);
	int n = get_n_eoprec(spacepos, timepos);
	int clerr = CL_SUCCESS;
	int evenodd = glob_pos % 2;

	(*latimer).reset();
	//CP: this is different than the host code, where this is done implicitly when converting the normal source to even/odd
	if(evenodd == 1){
		clerr = clSetKernelArg(create_point_source_eoprec,0,sizeof(cl_mem),&clmem_source_odd);
		if(clerr!=CL_SUCCESS) {
		  cout<<"clSetKernelArg 0 failed with " << clerr << " , aborting..."<<endl;
		  exit(HMC_OCLERROR);
		}
  }
	else {
		clerr = clSetKernelArg(create_point_source_eoprec,0,sizeof(cl_mem),&clmem_tmp_eoprec_2);
		if(clerr!=CL_SUCCESS) {
		  cout<<"clSetKernelArg 0 failed with " << clerr << " , aborting..."<<endl;
    exit(HMC_OCLERROR);
		}		
	}
	clerr = clSetKernelArg(create_point_source_eoprec,1,sizeof(int),&i);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source_eoprec,2,sizeof(int),&n);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,create_point_source_eoprec,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue create_point_source_eoprec kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
 	clFinish(queue);
	(*latimer).add();
	
	M_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_1, ls, gs, Mdiagtimer);
	dslash_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_3, EVEN, ls, gs, dslashtimer);

	saxpy_eoprec_device(clmem_tmp_eoprec_2, clmem_tmp_eoprec_3, clmem_one, clmem_source_even, ls, gs, latimer);
	
	return HMC_SUCCESS;
}


//functions to calculate the correlator
hmc_error Opencl_fermions::set_correlator_field_zero_device(const size_t ls, const size_t gs, usetimer * latimer){
  set_zero_spinorfield_device(clmem_corr, ls, gs, latimer);
  return HMC_SUCCESS;
}

hmc_error Opencl_fermions::add_solution_to_correlator_field_device(const size_t ls, const size_t gs, usetimer * latimer){
  saxpy_device(clmem_inout, clmem_corr, clmem_minusone, clmem_corr, ls, gs, latimer);
  return HMC_SUCCESS;
}

hmc_error Opencl_fermions::ps_correlator_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
	size_t ls = local_work_size;
	size_t gs = global_work_size;
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(ps_correlator,0,sizeof(cl_mem),&clmem_corr); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,ps_correlator,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue ps_correlator kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 clFinish(queue);
	
	(*timer).add();
	
	return HMC_SUCCESS;
}

cl_mem Opencl_fermions::get_clmem_inout(){
	return clmem_inout;
}

cl_mem Opencl_fermions::get_clmem_tmp(){
	return clmem_tmp;
}

cl_mem Opencl_fermions::get_clmem_inout_eoprec(){
	return clmem_inout_eoprec;
}

cl_mem Opencl_fermions::get_clmem_tmp_eoprec_1(){
	return clmem_tmp_eoprec_1;
}

hmc_error Opencl_fermions::finalize_fermions(){
	
  if(clReleaseKernel(ps_correlator)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(set_spinorfield_cold)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(M)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(saxpy)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(saxsbypz)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(scalar_product)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(scalar_product_reduction)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseKernel(set_zero_spinorfield)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(global_squarenorm)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(global_squarenorm_reduction)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(ratio)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(product)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(convert_to_kappa_format)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(convert_from_kappa_format)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(create_point_source)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(get_parameters()->get_use_eo()==TRUE) {
    if(clReleaseKernel(convert_to_kappa_format_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(convert_from_kappa_format_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(dslash_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(saxpy_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(scalar_product_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(set_zero_spinorfield_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(global_squarenorm_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseKernel(create_point_source_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
  	if(clReleaseKernel(M_sitediagonal)!=CL_SUCCESS) exit(HMC_OCLERROR);
  	if(clReleaseKernel(M_inverse_sitediagonal)!=CL_SUCCESS) exit(HMC_OCLERROR);
	}

  if(clReleaseMemObject(clmem_inout)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_source)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_rn)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_rhat)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_v)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_p)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_s)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_t)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_aux)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_tmp)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(get_parameters()->get_use_eo()==TRUE) {
    if(clReleaseMemObject(clmem_inout_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_source_even)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_source_odd)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_rn_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_rhat_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_v_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_p_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_s_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_t_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_aux_eoprec)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_tmp_eoprec_1)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_tmp_eoprec_2)!=CL_SUCCESS) exit(HMC_OCLERROR);
    if(clReleaseMemObject(clmem_tmp_eoprec_3)!=CL_SUCCESS) exit(HMC_OCLERROR);
  }  

  if(clReleaseMemObject(clmem_rho)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_rho_next)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_alpha)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_omega)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_beta)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_tmp1)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_tmp2)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_one)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_minusone)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_scalar_product_buf_glob)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_global_squarenorm_buf_glob)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_resid)!=CL_SUCCESS) exit(HMC_OCLERROR);
  /** @todo CP: this call seems to cause a seg fault!! */
  if(clReleaseMemObject(clmem_trueresid)!=CL_SUCCESS) exit(HMC_OCLERROR);

  return HMC_SUCCESS;
}

