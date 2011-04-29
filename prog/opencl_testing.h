/** @file
 * OpenCL test code
 */

#ifdef _FERMIONS_

//!!CP: this is here because the kernel global_squarenorm has the same name as the function in operations_spinor
// hmc_float global_squarenorm(hmc_spinor_field *field) {
hmc_float inline global_squarenorm_host(hmc_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<NTIME; t++) {
    for (int n=0; n<VOLSPACE; n++) {
      sum += local_squarenorm(field,n,t);
    }
  }
  return sum;
}
hmc_float inline global_squarenorm_eoprec_host(hmc_eoprec_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<EOPREC_SPINORFIELDSIZE; t++){
		
      sum += field[t].re*field[t].re + field[t].im*field[t].im;
  }
  return sum;
}

hmc_complex inline scalar_product_host (hmc_spinor_field * a, hmc_spinor_field * b){
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    res.re += a[n].re*b[n].re + a[n].im*b[n].im;
    res.im += a[n].re*b[n].im - a[n].im*b[n].re;
  	}
	return res;
}
	


#ifdef _TESTING_
hmc_error opencl::testing(hmc_gaugefield * gaugefield){
  cl_int clerr=CL_SUCCESS;

  cout<<"Begin testing OpenCL functions"<<endl;
  
  cout<<"Create test kernel..."<<endl;
  cl_kernel testkernel = clCreateKernel(clprogram,"test",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  int nsteps = 10;
  hmc_float beta = 4.2;

  //CP: this is made for a 4^4 lattice, otherwise one will run into problems with the VOL4D/2 definition!!!
  const size_t local_work_size  = 0;
  const size_t global_work_size = VOL4D/2;
  const size_t * local_work_size_p = (local_work_size == 0) ? 0 : &local_work_size;
 
  //CP: random number test
  hmc_float check=1;
  int size_1 = 3*1000;
  int size_2 = 3000;
  
  clmem_random_field_int = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(int)*size_1,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_random_field_float = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(float)*size_2,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_random_field_su2 = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float)*4,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_check = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_check,CL_TRUE,0,sizeof(hmc_float),&check,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  
  //CP: heatbath test (no reconstruct_twelve!!): set input matrices taken from a host_configuration
  hmc_su3matrix host_heatbath_link;
  (host_heatbath_link[0][0]).re = -0.349937;
  (host_heatbath_link[0][0]).im = -0.474231;
  (host_heatbath_link[0][1]).re = 0.301203;
  (host_heatbath_link[0][1]).im = -0.515532;
  (host_heatbath_link[0][2]).re = 0.544024;
  (host_heatbath_link[0][2]).im = -0.013786;
  (host_heatbath_link[1][0]).re = -0.789803;
  (host_heatbath_link[1][0]).im = -0.009681;
  (host_heatbath_link[1][1]).re = 0.047851;
  (host_heatbath_link[1][1]).im = 0.474619;
  (host_heatbath_link[1][2]).re = -0.083667;
  (host_heatbath_link[1][2]).im = 0.376251;  
  (host_heatbath_link[2][0]).re = 0.136194;
  (host_heatbath_link[2][0]).im = 0.101084 ;  
  (host_heatbath_link[2][1]).re = -0.637513;
  (host_heatbath_link[2][1]).im = -0.097608;  
  (host_heatbath_link[2][2]).re = 0.451216;
  (host_heatbath_link[2][2]).im = 0.593032;  

  hmc_staplematrix host_heatbath_staple;
  (host_heatbath_staple[0][0]).re = -0.043457;
  (host_heatbath_staple[0][0]).im = 2.350362;
  (host_heatbath_staple[0][1]).re = -2.827883;
  (host_heatbath_staple[0][1]).im = -0.147794;
  (host_heatbath_staple[0][2]).re = 0.239220;
  (host_heatbath_staple[0][2]).im = -0.353712;
  (host_heatbath_staple[1][0]).re = 2.094735;
  (host_heatbath_staple[1][0]).im = 2.500042;
  (host_heatbath_staple[1][1]).re = 0.274546;
  (host_heatbath_staple[1][1]).im = -2.518978;
  (host_heatbath_staple[1][2]).re = -1.411352;
  (host_heatbath_staple[1][2]).im = -0.568307;  
  (host_heatbath_staple[2][0]).re = 1.217603;
  (host_heatbath_staple[2][0]).im = 0.463087;  
  (host_heatbath_staple[2][1]).re = -0.337425;
  (host_heatbath_staple[2][1]).im = -0.909408;  
  (host_heatbath_staple[2][2]).re = 2.763673;
  (host_heatbath_staple[2][2]).im = -2.460356; 
  
  //set up rnd-array 
  int heatbath_rnd_array_size = 1e4;
  hmc_float * heatbath_rnd_array = (hmc_float*) malloc(sizeof(hmc_float)*heatbath_rnd_array_size);
  
  Random ran(50000);
  for (int i = 0; i<heatbath_rnd_array_size; i++){
    heatbath_rnd_array[i] = ran.doub();    
  }
  
  hmc_su3matrix host_heatbath_out;
  int heatbath_host_cter;
  //CP: there are 3 different possiblities to test the heatbath, with more or less random numbers
//   testing_heatbath_norandommat_no123(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2);
//   testing_heatbath_no123(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2, heatbath_rnd_array, &heatbath_host_cter);
   testing_heatbath(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2, heatbath_rnd_array, &heatbath_host_cter);
   
  clmem_heatbath_test_link_in = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_staple_in = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_link_out = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_rnd_array = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float)*heatbath_rnd_array_size, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_cter = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  hmc_complex tmp[9], tmp2[9];
  int i, j;
  for(i=0; i<NC; i++){
    for(j=0; j<NC; j++){
      tmp[j*NC+i].re = host_heatbath_link[i][j].re;
      tmp[j*NC+i].im = host_heatbath_link[i][j].im;
      tmp2[j*NC+i].re = host_heatbath_staple[i][j].re;
      tmp2[j*NC+i].im = host_heatbath_staple[i][j].im;     
      
  }}
  
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_link_in,CL_TRUE,0,sizeof(hmc_complex)*9,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_staple_in,CL_TRUE,0,sizeof(hmc_complex)*9,&tmp2,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_rnd_array,CL_TRUE,0,sizeof(hmc_float)*heatbath_rnd_array_size,heatbath_rnd_array,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
    
  //CP: solver test
  
  hmc_ocl_gaugefield* gfout = (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));
  copy_to_ocl_format(gfout, gaugefield);
  cl_mem clmem_gfout = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_gaugefield),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
  clerr = clEnqueueWriteBuffer(queue,clmem_gfout,CL_TRUE,0,sizeof(hmc_gaugefield),gfout,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
 
  //prepare simple fermionfield
  hmc_spinor_field solver_test_in[SPINORFIELDSIZE];
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(solver_test_in, n, t, a, j);
	}
      }
    }
  }
  hmc_float norm = global_squarenorm(solver_test_in);
  norm = sqrt(norm);
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    solver_test_in[n].re /= norm;
    solver_test_in[n].im /=norm;
  }
  
  int test_correlator_size = sizeof(hmc_complex)*NSPACE;
  int test_spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;

  cl_mem clmem_solver_test_spinor_in = clCreateBuffer(context,CL_MEM_READ_WRITE,test_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_solver_test_spinor_out = clCreateBuffer(context,CL_MEM_READ_WRITE,test_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_solver_test_correlator = clCreateBuffer(context,CL_MEM_READ_WRITE,test_correlator_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

  clerr = clEnqueueWriteBuffer(queue,clmem_solver_test_spinor_in,CL_TRUE,0,test_spinorfield_size,solver_test_in,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
  
  //CP: set kernel arguments
  
  clerr = clSetKernelArg(testkernel,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,2,sizeof(int),&nsteps);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,3,sizeof(cl_mem),&clmem_check);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,4,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,5,sizeof(cl_mem),&clmem_random_field_int);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,6,sizeof(cl_mem),&clmem_random_field_float);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,7,sizeof(cl_mem),&clmem_random_field_su2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 7 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,8,sizeof(int),&size_1);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 8 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,9,sizeof(int),&size_2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 9 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,10,sizeof(cl_mem),&clmem_heatbath_test_link_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 10 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,11,sizeof(cl_mem),&clmem_heatbath_test_staple_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 11 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,12,sizeof(cl_mem),&clmem_heatbath_test_link_out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 12 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,13,sizeof(cl_mem),&clmem_heatbath_test_rnd_array);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 13 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,14,sizeof(cl_mem),&clmem_heatbath_test_cter);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 14 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,15,sizeof(cl_mem),&clmem_gfout);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 15 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,16,sizeof(cl_mem),&clmem_solver_test_spinor_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 16 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,17,sizeof(cl_mem),&clmem_solver_test_spinor_out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 17 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,18,sizeof(cl_mem),&clmem_solver_test_correlator);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 18 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  cout << "Enqueue kernel.." << endl;
  clerr = clEnqueueNDRangeKernel(queue,testkernel,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 
  //CP: random number test
  
  int * random_int = (int*) malloc(sizeof(int)*size_1);
  float * random_float = (float*) malloc(sizeof(float)*size_2);
  hmc_float su2mat[4];

  clerr = clEnqueueReadBuffer(queue,clmem_random_field_int,CL_TRUE,0,sizeof(int)*size_1,random_int,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_random_field_float,CL_TRUE,0,sizeof(float)*size_2,random_float,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_random_field_su2,CL_TRUE,0,sizeof(hmc_float)*4,su2mat,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  
  //print result of test from LZ
  cout<<"\ttest functions: result: "<<check<<endl;
  
  //CP: print results from random number test
  cout << "\ttest random 1,2,3 (2 - 1/3*average of "<<size_1/3 << " per array position)" <<endl;
  float rnd_check_int_1 = 0;
  float rnd_check_int_2 = 0;
  float rnd_check_int_3 = 0;
  for(int i=0;i<size_1/3;i++){
// 	  cout << random_int[i*3] << "  " << random_int[i*3 + 1] << "  " << random_int[i*3 + 2] << " " ;
	rnd_check_int_1 += random_int[i*3 + 0];
	rnd_check_int_2 += random_int[i*3 + 1];
	rnd_check_int_3 += random_int[i*3 + 2];
  }
  rnd_check_int_1 = rnd_check_int_1/(size_1/3);
  rnd_check_int_2 = rnd_check_int_2/(size_1/3);
  rnd_check_int_3 = rnd_check_int_3/(size_1/3);
  cout << "\tcheck: " << 2. - rnd_check_int_1 << "  " << 2. - rnd_check_int_2 << "  " << 2. - rnd_check_int_3 <<  endl;
  cout << "\ttest random floats: (0.5 - average of "<<size_2 << ")" << endl;
  float rnd_check_float = 0;
  for(int i=0;i<size_2;i++){
    rnd_check_float += random_float[i];
  }
  
  rnd_check_float/=size_2;
  cout << "\tcheck: " << rnd_check_float-0.5 << endl;
  cout << "\ttest su2mat:";
  cout << "\t" << su2mat[0] << "  " << su2mat[1] << "  " << su2mat[2] << "  " << su2mat[3] << endl;
  cout << "\t1 - det is: " << 1. - (su2mat[0]*su2mat[0]+ su2mat[1]*su2mat[1] + su2mat[2]*su2mat[2] + su2mat[3]*su2mat[3]) << endl;

  
  //CP: heatbath_test
  clerr = clEnqueueReadBuffer(queue,clmem_heatbath_test_link_out,CL_TRUE,0,sizeof(hmc_complex)*9,tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  int heatbath_device_cter;
  clerr = clEnqueueReadBuffer(queue,clmem_heatbath_test_cter,CL_TRUE,0,sizeof(int),&heatbath_device_cter,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout << "\tresult of heatbath test: "<< endl;
  for(i = 0; i<NC; i++){
    for(j = 0; j<NC; j++){
      cout << "\t(" << (host_heatbath_out[i][j]).re - (tmp[i + NC*j]).re << ","  << 
              (host_heatbath_out[i][j]).im  - (tmp[i + NC*j]).im << ")";
    }cout << endl;}cout << endl;
  cout << "\thost code needed " << heatbath_host_cter << " random numbers" << endl;
  cout << "\tdevice code needed " << heatbath_host_cter - heatbath_device_cter << " more random numbers than host code" << endl;
  
  //CP: solver test
  //print result 
  
  hmc_complex solver_test_correlator[NSPACE];
  cout << "\ton the host the correlator is:" << endl;

  hmc_gaugefield * host_gaugefield;
  host_gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
  set_gaugefield_cold(host_gaugefield);
  simple_correlator(host_gaugefield, 0.125, 0.06, 0., 0., 0.,1000);

  clerr = clEnqueueReadBuffer(queue,clmem_solver_test_correlator,CL_TRUE,0,test_correlator_size,&solver_test_correlator,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "\tthe correlator is: " << endl;
  for(int z=0; z<NSPACE; z++) {
    printf("\t%d\t(%e,%e)\n",z, solver_test_correlator[z].re, solver_test_correlator[z].im);
  }
  
  //CP: finish
  if(clReleaseMemObject(clmem_check)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_float)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_int)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_su2)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_link_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_staple_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_link_out)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_cter)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_gfout)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_spinor_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_spinor_out)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_correlator)!=CL_SUCCESS) exit(HMC_OCLERROR);

  free(gfout);
  
  printf("opencl testing finished\n");
  return HMC_SUCCESS;
}
#endif //_TESTING_


hmc_error opencl::testing_spinor(inputparameters* parameters, size_t local_size, size_t global_size){

	usetimer noop;
	hmc_float norm;
	cout << "testing spinor kernels..." << endl;

	
	//set up spinor field for testing                                                                                                                          
	hmc_spinor_field test[SPINORFIELDSIZE];
	hmc_spinor_field test2[SPINORFIELDSIZE];
	init_spinorfield_cold(test);

  cout << "\tspinorfield norm: "<< global_squarenorm_host(test) << endl;
  copy_spinorfield_to_device(test, &noop);
  get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after copying to device: "<< global_squarenorm_host(test) << endl;
  copy_spinor_device(clmem_inout, clmem_v, &noop);
	set_zero_spinorfield_device(clmem_inout, local_size, global_size, &noop);
	get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after set_zero_spinorfield: "<< global_squarenorm_host(test) << endl;
	copy_spinor_device(clmem_v, clmem_inout, &noop);
	get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after copying spinorfield on device: "<< global_squarenorm_host(test) << endl;
	copy_float_from_device(clmem_resid, &norm, &noop);
	cout<< "\tglobal squarenorm on the device before calculation is: " << norm << endl;
	set_float_to_global_squarenorm_device(clmem_inout, clmem_resid, local_size, global_size, &noop);
	copy_float_from_device(clmem_resid, &norm, &noop);
	cout<< "\tglobal squarenorm on the device is: " << norm << endl;
	hmc_complex tester;
	copy_complex_from_device(clmem_tmp1, &tester, &noop);
	cout<< "\ttmp1 before copying is: " << tester.re << "," << tester.im << endl;	
	copy_complex_device(clmem_minusone, clmem_tmp1, &noop);
	copy_complex_from_device(clmem_tmp1, &tester, &noop);
	cout<< "\ttmp1 after copying minusone on it is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_ratio_device(clmem_tmp1, clmem_one, clmem_tmp2, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\ttmp2 after copying is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_product_device(clmem_tmp1, clmem_minusone, clmem_tmp2, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\ttmp2 after copying is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_scalar_product_device(clmem_inout, clmem_v, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(inout,inout) is: " << tester.re << "," << tester.im << endl;	
	
	saxpy_device(clmem_inout, clmem_v, clmem_minusone, clmem_rhat, local_size, global_size, &noop);
  set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t|(inout + inout)| is: " << tester.re << "," << tester.im << endl;
	
	saxsbypz_device(clmem_inout, clmem_v, clmem_inout, clmem_minusone, clmem_one, clmem_rhat, local_size, global_size, &noop);
  set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t|(inout + inout - inout)| is: " << tester.re << "," << tester.im << endl;
	
	cout << "\tcompare M on host and device..." << endl;
	
	hmc_gaugefield gaugefield;
	set_gaugefield_cold(&gaugefield);
	copy_gaugefield_to_device(&gaugefield, &noop);
	
	hmc_float kappa = (*parameters).get_kappa();
	hmc_float mu = (*parameters).get_mu();
	
	cout << "\tM on device:" ;
	
	init_spinorfield_cold(test);
	copy_spinorfield_to_device(test, &noop);
	
	M_device(clmem_inout, clmem_v, local_size, global_size, &noop, &noop, &noop);
	
	set_complex_to_scalar_product_device(clmem_v, clmem_v, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(Minout, Minout): " << tester.re << "," << tester.im << endl;	

	cout << "\tM on host:";
	
	M(test,test2,&gaugefield,kappa,mu, 0., 0., 0.);
	
	cout << "\t(M test2, M test2): "<< global_squarenorm_host(test2) << endl;
	
	cout << "\tcompare M_diag on host and device.." << endl;
	
	cout << "\tM_diag on dev:";
	
	init_spinorfield_cold(test);
	copy_spinorfield_to_device(test, &noop);

// 	copied from M_device
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(M_diag,0,sizeof(cl_mem),&clmem_inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,1,sizeof(cl_mem),&clmem_inout);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,2,sizeof(cl_mem),&clmem_kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,3,sizeof(cl_mem),&clmem_mu);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M_diag,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue M_diag kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	set_complex_to_scalar_product_device(clmem_inout, clmem_inout, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(M_diag inout, M_diag inout): " << tester.re << "," << tester.im << endl;	
	
	cout << "\tM_diag on host:";
// 		M_diag(test, test2, kappa, mu);  
// 	copied from host_operations_fermionmatrix since M_diag is also a kernel
	hmc_spinor spinout[SPINORSIZE];
	//iterate over all lattice sites
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			get_spinor_from_field(test,spinout,spacepos,timepos);
			M_diag_local(spinout, kappa, mu);
			put_spinor_to_field(spinout,test2,spacepos,timepos);
		}}
		
	cout << "\t(M_diag test, M_diag test): "<< global_squarenorm_host(test2) << endl;
	
		set_gaugefield_cold(&gaugefield);
	copy_gaugefield_to_device(&gaugefield, &noop);
	
	cout << "\tdslash on dev:";
	
	init_spinorfield_cold(test);
	copy_spinorfield_to_device(test, &noop);
	
	clerr = clSetKernelArg(dslash,0,sizeof(cl_mem),&clmem_inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,1,sizeof(cl_mem),&clmem_v);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,2,sizeof(cl_mem),&clmem_gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,3,sizeof(cl_mem),&clmem_theta_fermion);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,4,sizeof(cl_mem),&clmem_chem_pot_re);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,5,sizeof(cl_mem),&clmem_chem_pot_im);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,dslash,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	
	set_complex_to_scalar_product_device(clmem_v, clmem_v, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);

	cout<< "\t(dslash inout, dslash inout): " << tester.re << "," << tester.im << endl;	
	
	cout << "\tdslash on host:";
	init_spinorfield_cold(test);
// 		dslash(test,test,gaugefield, 0., 0., 0.,);
// 	copied from host_operations_fermionmatrix since M_diag is also a kernel
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			set_local_zero_spinor(spinout);   
			//like in host_geometry
			int coord[NDIM];
			coord[0]=0;
			for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
			
			// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
			dslash_temporal(spinout, spacepos, timepos, test, &gaugefield, 0., 0., 0.);
			// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
			dslash_spatial (spinout, coord, 1, spacepos, timepos, test, &gaugefield, 0., 0., 0.);
			// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
			dslash_spatial (spinout, coord, 2, spacepos, timepos, test, &gaugefield, 0., 0., 0.);
			// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
			dslash_spatial (spinout, coord, 3, spacepos, timepos, test, &gaugefield, 0., 0., 0.);
			
			put_spinor_to_field(spinout,test2,spacepos,timepos);
    }
  }
  
	tester = scalar_product_host(test2, test2);
	cout << "\t(dslash test, dslash test): "<< tester.re << " , " << tester.im << endl;
	
	
	hmc_spinor_field b[SPINORFIELDSIZE];
	//CP: a gaugefield is not needed here
	hmc_gaugefield dummy;
	set_gaugefield_cold(&dummy);
// 	create_point_source(b,1,0,0,0.15,4.,&dummy);
	//CP: taken from create_point_source
  for (int i = 0; i< SPINORFIELDSIZE; i++) {b[i].re = 0.; b[i].im = 0.;}

  int color = spinor_color(1);
  int spin = spinor_spin(1,color);

  b[spinor_field_element(spin,color,0,0)].re = sqrt(2.*kappa);

	copy_source_to_device(b, &noop);
	set_complex_to_scalar_product_device(clmem_source, clmem_source, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source, source) is: " << tester.re << "," << tester.im << endl;	
	
	cout << "\tperform bicgstab..." << endl;
	bicgstab_device(&noop, &noop, &noop, &noop, &noop, &noop,&noop,  local_size, global_size, 100);
	
	
	cout<< "calculate simple_correlator with device-solver..." << endl;
	hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }

	init_spinorfield_cold(test);
	copy_spinorfield_to_device(test, &noop);


	int cgmax = 100;
	hmc_spinor_field phi[SPINORFIELDSIZE];
  for(int k=0; k<NC*NSPIN; k++) {
		hmc_spinor_field b[SPINORFIELDSIZE];
		//create_point_source(b,k,0,0,kappa,mu,&gaugefield);
		//CP: taken from create_point_source
		for (int i = 0; i< SPINORFIELDSIZE; i++) {b[i].re = 0.; b[i].im = 0.;}

		int color = spinor_color(k);
		int spin = spinor_spin(k,color);

		b[spinor_field_element(spin,color,0,0)].re = sqrt(2.*kappa);
	
		copy_source_to_device(b, &noop);
		convert_to_kappa_format_device(clmem_inout, local_size, global_size, &noop);
		bicgstab_device(&noop, &noop, &noop, &noop, &noop,&noop,&noop,  local_size, global_size, cgmax);
		convert_from_kappa_format_device(clmem_inout, clmem_inout, local_size, global_size, &noop);
		get_spinorfield_from_device(phi, &noop);

		for(int timepos = 0; timepos<NTIME; timepos++) {
			for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
				for(int alpha = 0; alpha<NSPIN; alpha++) {
					for(int c = 0; c<NC; c++) {
					// int j = spinor_element(alpha,c);
					int n = spinor_field_element(alpha, c, spacepos, timepos);
					int z = get_spacecoord(spacepos, 3);
					hmc_complex tmp = phi[n];
					hmc_complex ctmp = complexconj(&tmp);
					hmc_complex incr = complexmult(&ctmp,&tmp);
					correlator_ps[z].re += incr.re;
					correlator_ps[z].im += incr.im;
		}}}}
  }

  printf("\tpseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }

	
	cout << "\ttesting eoprec functions..." << endl;

	cout << "\tset sources to zero:" << endl;
	set_zero_spinorfield_eoprec_device(clmem_source_even, local_size, global_size, &noop);
	set_zero_spinorfield_eoprec_device(clmem_source_odd, local_size, global_size, &noop);
	set_complex_to_scalar_product_eoprec_device(clmem_source_even, clmem_source_even, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_even, source_even) is: " << tester.re << "," << tester.im << endl;
	set_complex_to_scalar_product_eoprec_device(clmem_source_odd, clmem_source_odd, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_odd, source_odd) is: " << tester.re << "," << tester.im << endl;
	
	cout << "\tcreate pointsource:" << endl;
	create_point_source_eoprec_device(1,0,0,local_size, global_size, &noop, &noop, &noop);
	set_complex_to_scalar_product_eoprec_device(clmem_source_even, clmem_source_even, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_even, source_even) is: " << tester.re << "," << tester.im << endl;
	set_complex_to_scalar_product_eoprec_device(clmem_source_odd, clmem_source_odd, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_odd, source_odd) is: " << tester.re << "," << tester.im << endl;
	
	cout << "\tperform saxpy: " << endl;
	saxpy_eoprec_device(clmem_source_odd, clmem_source_even, clmem_one, clmem_source_even, local_size, global_size, &noop);
	set_complex_to_scalar_product_eoprec_device(clmem_source_even, clmem_source_even, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_even, source_even) is: " << tester.re << "," << tester.im << endl;
	
	cout << "\tperform saxsbypz: " << endl;
	saxsbypz_eoprec_device(clmem_source_odd, clmem_source_even, clmem_source_odd, clmem_kappa, clmem_one, clmem_source_even, local_size, global_size, &noop);
	set_complex_to_scalar_product_eoprec_device(clmem_source_even, clmem_source_even, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source_even, source_even) is: " << tester.re << "," << tester.im << endl;
	
	cout << "\tperform global_squarenorm" << endl;
	hmc_float resid;
	set_float_to_global_squarenorm_eoprec_device(clmem_source_even, clmem_resid, local_work_size, global_work_size, &noop);
	copy_float_from_device(clmem_resid, &resid, &noop);
	cout << "\t|source_even|^2: " << resid << endl;

		hmc_eoprec_spinor_field in_even [EOPREC_SPINORFIELDSIZE];
		hmc_eoprec_spinor_field in_odd [EOPREC_SPINORFIELDSIZE];
		
	init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);
	
	M_sitediagonal_device(clmem_inout_eoprec, clmem_source_odd, local_size, global_size, &noop);
	set_float_to_global_squarenorm_eoprec_device(clmem_source_odd, clmem_resid, local_work_size, global_work_size, &noop);
	copy_float_from_device(clmem_resid, &resid, &noop);
	cout << "\t|M_sitediag inout_even|^2 dev: " << resid << endl;	
	
		init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);
	
	//CP: from: hmc_error M_sitediagonal(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_float kappa, hmc_float mu){
	//iterate over half the lattice
	set_local_zero_spinor(spinout);
  for(int n=0; n<VOL4D/2; n++) {
    get_spinor_from_eoprec_field(in_even,spinout,n);
    M_diag_local(spinout, kappa, mu);
    put_spinor_to_eoprec_field(spinout,in_odd,n);
  }
	cout << "\t|M_sitediag inout_even|^2 host: " << global_squarenorm_eoprec_host(in_odd) << endl;	

	init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);
	
	M_inverse_sitediagonal_device(clmem_inout_eoprec, clmem_source_odd, local_size, global_size, &noop);
	set_float_to_global_squarenorm_eoprec_device(clmem_source_odd, clmem_resid, local_work_size, global_work_size, &noop);
	copy_float_from_device(clmem_resid, &resid, &noop);
	cout << "\t|M_inv_sitediag inout_even|^2 dev: " << resid << endl;	
	
	init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);
		
	//CP: from: hmc_error M_inverse_sitediagonal(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_float kappa, hmc_float mu){
	set_local_zero_spinor(spinout);
	//iterate over half the lattice
	for(int n=0; n<VOL4D/2; n++) {
		hmc_float minuskappa = -kappa;
		get_spinor_from_eoprec_field(in_even,spinout,n);
		M_diag_local(spinout, minuskappa, mu);
		hmc_float denom = 1. + 4.*kappa*kappa*mu*mu;
		real_multiply_spinor(spinout,1./denom);
		put_spinor_to_eoprec_field(spinout,in_odd,n);
	}
	
	cout << "\t|M_inv_sitediag inout_even|^2 host: " << global_squarenorm_eoprec_host(in_odd) << endl;	
	
	dslash_eoprec_device(clmem_inout_eoprec, clmem_source_odd, ODD, local_size, global_size, &noop);
	set_float_to_global_squarenorm_eoprec_device(clmem_source_odd, clmem_resid, local_work_size, global_work_size, &noop);
	copy_float_from_device(clmem_resid, &resid, &noop);
	cout << "\t|dslash inout_even|^2 dev: " << resid << endl;	
	
	init_spinorfield_cold_eoprec(in_even);
	
	//CP: from hmc_error dslash_eoprec(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int evenodd){
	set_local_zero_spinor(spinout);
	int evenodd = ODD;
	int ns, nt;
	//iterate over half the lattice
  for(int n=0; n<VOL4D/2; n++) {
    set_local_zero_spinor(spinout);    

    if(evenodd == ODD) get_odd_site(n, &ns, &nt);
    else get_even_site(n, &ns, &nt);
       
    //like in host_geometry
    int coord[NDIM];
    coord[0]=0;
    for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(ns,j);    
    
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal_eoprec (spinout, ns, nt, in_even, &dummy, 0.,0.,0.);
		// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
		dslash_spatial_eoprec (spinout, coord, 1, ns, nt, in_even, &dummy, 0.,0.,0.);
		// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
		dslash_spatial_eoprec (spinout, coord, 2, ns, nt, in_even, &dummy, 0.,0.,0.);   
		// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
		dslash_spatial_eoprec (spinout, coord, 3, ns, nt, in_even, &dummy, 0.,0.,0.);
		
		//!!CP: why this kappa??
		real_multiply_spinor(spinout,-kappa);
		put_spinor_to_eoprec_field(spinout,in_odd,n);
  }
	cout << "\t|dslash inout_even|^2 host: " << global_squarenorm_eoprec_host(in_odd) << endl;	
	
	cout << "\tperform bicgstab_eoprec..." << endl;
	
	
	init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);
	set_float_to_global_squarenorm_eoprec_device(clmem_inout_eoprec, clmem_resid, local_work_size, global_work_size, &noop);
	copy_float_from_device(clmem_resid, &resid, &noop);
	cout << "\t|inout_eoprec|^2: " << resid << endl;
  create_point_source_eoprec_device(0,0,0,local_size, global_size, &noop, &noop, &noop);
	
	bicgstab_eoprec_device(&noop, &noop, &noop, &noop, &noop,&noop, &noop, local_size, global_size, 1);

	cout<< "calculate simple_correlator with eoprec device-solver..." << endl;
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }

	init_spinorfield_cold_eoprec(in_even);
	copy_eoprec_spinorfield_to_device(in_even, &noop);

	cgmax = 100;
	hmc_eoprec_spinor_field phi_even[EOPREC_SPINORFIELDSIZE];
	hmc_eoprec_spinor_field phi_odd[EOPREC_SPINORFIELDSIZE];
  for(int k=0; k<NC*NSPIN; k++) {
		create_point_source_eoprec_device(k,0,0,local_size, global_size, &noop, &noop, &noop);

		//CP: even solution
		convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, local_size, global_size, &noop);
		bicgstab_eoprec_device(&noop, &noop, &noop, &noop, &noop, &noop, &noop,local_size, global_size, cgmax);	
		convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, local_size, global_size, &noop);
		get_eoprec_spinorfield_from_device(phi_even, &noop);

		//P: odd solution
		convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, local_size, global_size, &noop);
		dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_3, ODD, local_size, global_size, &noop);
		M_inverse_sitediagonal_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1,local_size, global_size, &noop);
		M_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_2,local_size, global_size, &noop);
		saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, clmem_one, clmem_inout_eoprec, local_size, global_size, &noop);
		convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, local_size, global_size, &noop);
		get_eoprec_spinorfield_from_device(phi_odd, &noop);

		//CP: whole solution
		convert_from_eoprec(phi_even,phi_odd,phi);
	
		for(int timepos = 0; timepos<NTIME; timepos++) {
			for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
				for(int alpha = 0; alpha<NSPIN; alpha++) {
					for(int c = 0; c<NC; c++) {
					// int j = spinor_element(alpha,c);
					int n = spinor_field_element(alpha, c, spacepos, timepos);
					int z = get_spacecoord(spacepos, 3);
					hmc_complex tmp = phi[n];
					hmc_complex ctmp = complexconj(&tmp);
					hmc_complex incr = complexmult(&ctmp,&tmp);
					correlator_ps[z].re += incr.re;
					correlator_ps[z].im += incr.im;
		}}}}
  }

  printf("\tpseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }

	cout << "...testing fermion kernels done" << endl;
	
	return HMC_SUCCESS;
}

#endif /* _FERMIONS_ */
