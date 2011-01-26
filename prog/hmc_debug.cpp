#include "hmc.h"

int main(int argc, char* argv[]) {

  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);
  print_info(&parameters);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
  hmc_rndarray rndarray;

  init_gaugefield(gaugefield,&parameters,&inittime);

  //this needs optimization
  const size_t local_work_size  = VOL4D/2;
  const size_t global_work_size = local_work_size;
  //one should define a definite number of threads and use this here
  init_random_seeds(rnd, rndarray, VOL4D/2, &inittime);
  
     hmc_su3matrix A;
    (A[0][0]).re = 0.042082;
    (A[0][0]).im = -0.203080;
    (A[0][1]).re =-0.911624;
    (A[0][1]).im =-1.019873;
    (A[0][2]).re =-0.541189;
    (A[0][2]).im =0.773009;
    (A[1][0]).re =-0.584805;
    (A[1][0]).im =-0.370283;
    (A[1][1]).re =0.066303;
    (A[1][1]).im =0.222891;
    (A[1][2]).re =0.239045;
    (A[1][2]).im =0.477723;
    (A[2][0]).re =0.435155;
    (A[2][0]).im =0.627326;
    (A[2][1]).re =0.042681;
    (A[2][1]).im =-0.258224;
    (A[2][2]).re =-0.162545;
    (A[2][2]).im =0.456923;
    
//     for(int i = 0; i<NC; i++){
// 	for(int j = 0; j<NC; j++){
// 	  A[i][j] = hmc_complex_zero;
//       }}
      
    cout << "adjoint A: " << endl;
    hmc_su3matrix dummy2, dummy3;
    copy_su3matrix(&dummy2, &A);
    copy_su3matrix(&dummy3, &A);
    adjoin_su3matrix(&dummy3);
     print_su3mat(&dummy2);
     print_su3mat(&dummy3);

    accumulate_su3matrix_prod(&dummy3, &dummy2);
//     adjoin_su3matrix(&dummy3);
    for(int i = 0; i<NC; i++){
	for(int j = 0; j<NC; j++){
	  cout << "(" << (dummy3[i][j]).re << "," << (dummy3[i][j]).im << ") ";
	}cout << endl;}cout << endl;
    
	    copy_su3matrix(&dummy2, &A);
    copy_su3matrix(&dummy3, &A);
    adjoin_su3matrix(&dummy3);
    accumulate_su3matrix_prod(&dummy2, &dummy3);
    for(int i = 0; i<NC; i++){
	for(int j = 0; j<NC; j++){
	  cout << "(" << (dummy2[i][j]).re << "," << (dummy2[i][j]).im << ") ";
	}cout << endl;}cout << endl;
	
    
    
    int pos=0, t=0, mu =0;
     put_su3matrix(gaugefield, &A, pos, t, mu);
//     mu = 2;
    pos=1;
//      put_su3matrix(gaugefield, &A, pos, t, mu);
  mu = 0;
  pos=0;
  
  
  opencl gpu(CL_DEVICE_TYPE_GPU, &inittime);

  cout << "initial values of observables:\n\t" ;
  print_gaugeobservables(gaugefield, &polytime, &plaqtime);

  cout << "make one update on the host..."<< endl;
//   heatbath_update (gaugefield, parameters.get_beta());
  cout << "now values of observables:\n\t" ;
  print_gaugeobservables(gaugefield, &polytime, &plaqtime);
  
  gpu.copy_gaugefield_to_device(gaugefield, &copytime);
  gpu.copy_rndarray_to_device(rndarray, &copytime);

  
  heatbath_overrelax(gaugefield, parameters.get_beta());
  hmc_su3matrix dummy;
  mu = 0;

  print_gaugeobservables(gaugefield, &polytime, &plaqtime);

  cout << "pos=0" << endl;
  pos=0;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  
  cout << "pos=1" << endl;
    pos=1;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  
  
//    gpu.testing();
  
  //this has go into a function later
  int nsteps = 1, writeout=nsteps, measure=int(nsteps/2);
  
  cout<<"perform "<<nsteps<<" heatbath steps on OpenCL device..."<<endl;
  for(int i = 0; i<nsteps; i++){
   // gpu.run_heatbath(parameters.get_beta(), local_work_size, global_work_size, &updatetime);
     gpu.run_overrelax(parameters.get_beta(), local_work_size, global_work_size, &overrelaxtime);
    //this has to be done on the GPU
//     if(((i+1)%writeout)== 0) gpu.get_gaugefield_from_device(gaugefield, &copytime);
//     if(((i+1)%measure) == 0) print_gaugeobservables(gaugefield, &polytime, &plaqtime);
  }

  gpu.get_gaugefield_from_device(gaugefield, &copytime);

  print_gaugeobservables(gaugefield, &polytime, &plaqtime);
 
  cout << "pos=0" << endl;
  pos=0;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  
  cout << "pos=1" << endl;
    pos=1;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy); 
  
  hmc_su3matrix U; hmc_complex det; hmc_complex detadj;hmc_complex detsqnorm;
  
  cout << " t=1" << endl;
  t=1;
    cout << "pos=0" << endl;
  pos=0;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  
  cout << "pos=1" << endl;
    pos=1;
  mu = 0;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 1;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 2;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy);  
  mu = 3;
  get_su3matrix(&dummy, gaugefield, pos, t, mu);
  print_su3mat(&dummy); 

  
  
  //test1
 /*  
  int pos = 0, t = 0, mu=3;
  get_su3matrix(&U, gaugefield, pos, t, mu);
  print_su3mat(&U);
  det = det_su3matrix(&U);
  cout << det.re << "  "  << det.im << endl;
  detadj = complexconj(&det);
//   cout << detadj.re << "  " << detadj.im << endl;
  detsqnorm = complexmult(&det, &detadj);
//   cout << detsqnorm.re << "  " << detsqnorm.im << endl;
  project_su3(&U);
  print_su3mat(&U);
  det = det_su3matrix(&U);
  cout << det.re << "  " << det.im << endl;

  cout << "project whole lattice: " << endl;
  for(int i = 0; i<2; i++){
  for(pos = 0; pos<VOLSPACE; pos++){
    for(t = 0; t<NTIME; t++){
      for(mu=0; mu<NDIM; mu++){
        get_su3matrix(&U, gaugefield, pos, t, mu);
        project_su3(&U);
        put_su3matrix(gaugefield, &U, pos, t, mu);
      }}}}
  print_gaugeobservables(gaugefield, &polytime, &plaqtime);

  
  cout << "look at one polyakovloop:" << endl;
  int const tdir = 0;
  hmc_complex res;
  res.re=0;
  res.im=0;
  int n = 0;
  //  for(int n=0; n<VOLSPACE; n++) {                                                                                                                        
    hmc_su3matrix prod;
    unit_su3matrix(&prod);
    for(int t=0; t<NTIME; t++) {
      hmc_su3matrix tmp;
      get_su3matrix(&tmp,gaugefield,n,t,tdir);
      print_su3mat(&tmp);
      accumulate_su3matrix_prod(&prod,&tmp);
      project_su3(&prod);
    }
    hmc_complex tmpcomplex = trace_su3matrix(&prod);
    complexaccumulate(&res,&tmpcomplex);
    //}                                         
    res.re /= static_cast<hmc_float>(NC);
    res.im /= static_cast<hmc_float>(NC);
    print_su3mat(&prod);
    det = det_su3matrix(&prod);
    cout << det.re << "  " << det.im << endl;
    project_su3(&prod);
    print_su3mat(&prod);
    det = det_su3matrix(&prod);
    cout << det.re << "  " << det.im << endl;

    cout << res.re << "  " << res.im << endl;
    */
 
 //Test2

/*
    cout << "Test overrelaxing.." << endl;

    
    print_su3mat(&A);
    
    cout << "make overrelaxation.." <<endl;

      hmc_staplematrix W,  staplematrix;
      hmc_complex w [su2_entries];
      hmc_float w_pauli[su2_entries], k;
      int order[3]; 
      //random_1_2_3(order);
      order[0] = 1; order[1] = 2; order[2] = 3;
      
//       det = det_su3matrix(&A);
//       detadj = complexconj(&det);
//       detsqnorm = complexmult(&det, &detadj);
	project_su3(&A);
      
      //calc_staple(gaugefield, &staplematrix, pos, t, mu);
      for(int i = 0; i<NC; i++){
	for(int j = 0; j<NC; j++){
	  if(i==j) {(staplematrix[i][j]).re = 6.;
	  	(staplematrix[i][j]).im = 0.;}
	  else staplematrix[i][j] = hmc_complex_zero;
      }}

      for(int i=0; i<1; i++)
      {
        //Produkt aus U und staple, saved in W
        multiply_staplematrix(&W, &A, &staplematrix);
// 	print_staplemat(&W);
	//Reduktion auf SU(2) submatrix
        reduction (w, W, order[i]);
	
        //Darstellung von W in Paulibasis
        w_pauli[0] = 0.5*(w[0].re + w[3].re);
        w_pauli[1] = 0.5*(w[1].im + w[2].im);
        w_pauli[2] = 0.5*(w[1].re - w[2].re);
        w_pauli[3] = 0.5*(w[0].im - w[3].im);
     
        //Berechnung der Normierung k, entspricht Determinante in normaler Basis
        k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
        hmc_float a[su2_entries];
//         a[0] = hmc_float(1)/k * w_pauli[0];
//         a[1] = hmc_float(-1)/k * w_pauli[1];
//         a[2] = hmc_float(-1)/k * w_pauli[2];
//         a[3] = hmc_float(-1)/k * w_pauli[3];
        a[0] = w_pauli[0];
        a[1] = -w_pauli[1];
        a[2] = -w_pauli[2];
        a[3] = -w_pauli[3];

        //Square a and save in w
        w_pauli[0] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
        w_pauli[1] = 2.*a[0]*a[1];
        w_pauli[2] = 2.*a[0]*a[2];
        w_pauli[3] = 2.*a[0]*a[3];
      
        //go back to a su2 matrix in standard basis
        w[0].re = w_pauli[0];
        w[0].im = w_pauli[3];
        w[1].re = w_pauli[2];
        w[1].im = w_pauli[1];
        w[2].re = -w_pauli[2];
        w[2].im = w_pauli[1];
        w[3].re = w_pauli[0];
        w[3].im = -w_pauli[3];
	
        //extend to SU3
        hmc_su3matrix extW, tmp;
        extend (&extW, order[i], w);
	
	print_staplemat(&extW);
	
	multiply_su3matrices(&tmp, &extW, &A);
	copy_su3matrix(&A, &tmp);
      }
//       put_su3matrix(gaugefield, &U, pos, t, mu);
 
      print_su3mat(&A);
    */
  
/*
* The whole non-opencl part

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Thermalization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(!parameters.get_startcondition()==START_FROM_SOURCE){
    cout << endl << "perform thermalization" << endl;
    for (int i = 0; i<parameters.get_thermalizationsteps(); i++){
      updatetime.reset();
      heatbath_update (gaugefield, parameters.get_beta());
      updatetime.add();
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Test-Measurements
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  printf("Initial measurements...\n");
  
  plaqtime.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  plaqtime.add();
  
  polytime.reset();
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  polytime.add();
  
  //Testing
  for (int i = 0; i<parameters.get_heatbathsteps()/2-1; i++){
    plaqtime.reset();
    plaquette(gaugefield);
    plaqtime.add();
    
    polytime.reset();
    polyakov(gaugefield);
    polyakov_x(gaugefield);
    polyakov_y(gaugefield);
    polyakov_z(gaugefield);
    polytime.add();
  }
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Gaugefield-updates in thermalized system on CPU
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  printf("\nheatbath...\n");
  
  
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
    updatetime.reset();
    heatbath_update (gaugefield, parameters.get_beta());
    updatetime.add();
  }
  
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
    overrelaxtime.reset();
    heatbath_overrelax (gaugefield, parameters.get_beta());
    overrelaxtime.add();
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // More measurements
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  printf("now the stuff is..\n");  
  
  plaqtime.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  plaqtime.add();
  
  polytime.reset();
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  polytime.add();
  
*/  
  
  totaltime.add();
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  time_output(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  free(gaugefield);

  return HMC_SUCCESS;
}
