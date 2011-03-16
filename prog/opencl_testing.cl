//opencl_testing.cl

void testing_heatbath_norandommat_no123(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //SU2Update(r_pauli, beta_neu);
    //random su2-matrix
    r_pauli[0] = 0.748102;
    r_pauli[1] = 0.293055;
    r_pauli[2] = 0.591048;
    r_pauli[3] = -0.0715786;

    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  
  return;
}

void testing_heatbath_no123(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  int cter = 0;
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;
    
    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);


    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

void testing_heatbath(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  int cter = 0;
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  
  order[0] = convert_int( rnd_array[cter] * 3 ) + 1;
  cter ++;
  do
    {order[1] = convert_int( rnd_array[cter] * 3 ) + 1;
     cter ++;}
  while (order[1] == order[0]);
  order[2] = 6 - order[1] - order[0];
    
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;
    
    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);


    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

__kernel void test(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int nsteps,__global hmc_float* check,__global hmc_ocl_gaugefield* gaugefield2, __global hmc_ocl_ran * rnd, __global int * random_field_int, __global float *  random_field_float, __global hmc_float * su2mat, const int size_1, const int size_2, __global hmc_ocl_su3matrix * heatbath_link_in, __global hmc_ocl_su3matrix * heatbath_staple_in, __global hmc_ocl_su3matrix * heatbath_link_out, __global hmc_float * heatbath_rnd_array_in,__global int * heatbath_cter)
{
>>>>>>> code/master
  
  return;
}

void testing_heatbath_no123(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
//   random_1_2_3(order);
  order[0] = 1; order[1] = 2; order[2] = 3;
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  int cter = 0;
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;
    
    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);


    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

void testing_heatbath(hmc_ocl_su3matrix * in, hmc_ocl_staplematrix * staple_in, hmc_ocl_su3matrix * out, hmc_float beta, hmc_float * rnd_array, int * cter_out){
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
      
  hmc_ocl_su3matrix su3_in[SU3SIZE];
  copy_su3matrix(su3_in, in);
  int cter = 0;
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  
  order[0] = convert_int( rnd_array[cter] * 3 ) + 1;
  cter ++;
  do
    {order[1] = convert_int( rnd_array[cter] * 3 ) + 1;
     cter ++;}
  while (order[1] == order[0]);
  order[2] = 6 - order[1] - order[0];
    
  hmc_complex det = det_su3matrix(su3_in);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(su3_in);
  
  for(int i=0; i<3; i++)
  {
    //Produkt aus U und staple, saved in W
    multiply_staplematrix(W, su3_in, staple_in);
    //Reduktion auf SU(2) submatrix
    reduction (w, W, order[i]);
    //Darstellung von W in Paulibasis
    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);

    //Berechnung der Normierung k, entspricht Determinante in normaler Basis
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
    w_pauli[0] = w_pauli[0]/k;
    w_pauli[1] = -w_pauli[1]/k;
    w_pauli[2] = -w_pauli[2]/k;
    w_pauli[3] = -w_pauli[3]/k;
	
    //beta' = 2Beta/Nc*k
    hmc_float beta_neu =  2.*beta / convert_float(NC)*k;
      
    //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
    //function  SU2Update(r_pauli, beta_neu):
    hmc_float delta;
    hmc_float a0 ;
    hmc_float eta ;
    hmc_float alpha = beta_neu;
    
    do
    {
      hmc_float rnd1 = rnd_array[cter];
      hmc_float rnd2 = rnd_array[cter+1];
      hmc_float rnd3 = rnd_array[cter+2];
      hmc_float rnd4 = rnd_array[cter+3];
      delta = -log(rnd1)/alpha*pow(cos(2. * PI * rnd2), 2.) -log(rnd3)/alpha;
      a0 = 1.-delta;
      eta = rnd4;
      cter += 4;
    }while ( (1.-0.5*delta) < eta*eta);
    hmc_float rnd5 = rnd_array[cter];
    hmc_float rnd6 = rnd_array[cter+1];
    cter += 2;
    hmc_float phi = 2.*PI*rnd5;
    hmc_float theta = asin(2.*rnd6 - 1.);
    r_pauli[0] = a0;
    r_pauli[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
    r_pauli[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
    r_pauli[3] = sqrt(1.-a0 * a0)*sin(theta);


    //Multipliziere neuen Link mit w, alles in Paulibasis
    hmc_float su2_tmp[su2_entries];
    su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
    su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
    su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
    su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
    r_pauli[0] = su2_tmp[0];      
    r_pauli[1] = su2_tmp[1];
    r_pauli[2] = su2_tmp[2];
    r_pauli[3] = su2_tmp[3];

      
    //go back to a su2 matrix in standard basis
    w[0].re = r_pauli[0];
    w[0].im = r_pauli[3];
    w[1].re = r_pauli[2];
    w[1].im = r_pauli[1];
    w[2].re = -r_pauli[2];
    w[2].im = r_pauli[1];
    w[3].re = r_pauli[0];
    w[3].im = -r_pauli[3];
 
    //extend to SU3
    hmc_ocl_su3matrix extW[SU3SIZE], tmp[SU3SIZE];
    extend (extW, order[i], w);
	
    multiply_su3matrices(tmp, extW, su3_in);
    copy_su3matrix(su3_in, tmp);
  }
  copy_su3matrix(out, su3_in);
  (*cter_out) = cter;
  return;
}

__kernel void test(
__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int nsteps,__global hmc_float* check
//random_test-args
, __global hmc_ocl_ran * rnd, __global int * random_field_int, __global float * random_field_float, __global hmc_float * su2mat, const int size_1, const int size_2
//heatbath_test-args
 ,__global hmc_ocl_su3matrix * heatbath_link_in, __global hmc_ocl_su3matrix * heatbath_staple_in, __global hmc_ocl_su3matrix * heatbath_link_out, __global hmc_float * heatbath_rnd_array_in,__global int * heatbath_cter
//solver_test-args
,__global hmc_ocl_gaugefield* gaugefield2, __global hmc_spinor_field * solver_spinor_in, __global hmc_spinor_field * solver_spinor_out, __global hmc_complex * solver_correlator
)
{
  int id = get_global_id(0);
  if(id >0) return;
  else{
  //test by LZ
  hmc_complex testsum;
  testsum.re = 0;
  testsum.im = 0;
  hmc_complex ctmp;
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_staplematrix tester[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix tester2[STAPLEMATRIXSIZE];
  unit_su3matrix(prod);
  if(id==0) {
  for(int spacepos = 0; spacepos < VOLSPACE; spacepos++) {
    for(int mu = 0; mu<NDIM; mu++) {
      for(int t=0; t<NTIME; t++) {
	hmc_ocl_su3matrix tmp[SU3SIZE];
	get_su3matrix(tmp,gaugefield,spacepos,t,mu);
	hmc_ocl_su3matrix tmp2[SU3SIZE];
	copy_su3matrix(tmp2,tmp);

	//test matrix operations...
	hmc_ocl_su3matrix tmp3[SU3SIZE];
	multiply_su3matrices(tmp3,tmp2,tmp);
	adjoin_su3matrix(tmp3);
	accumulate_su3matrix_prod(prod, tmp3);
	ctmp = det_su3matrix(tmp3);
	ctmp.re -= hmc_one_f;
	hmc_complex ctmpconj = complexconj(&ctmp);
	hmc_complex square = complexmult(&ctmp,&ctmpconj);
	complexaccumulate(&testsum,&square);

        copy_staplematrix(tester, tester2);

	//	unit_su3matrix(tmp2);
	//	zero_su3matrix(tmp2);

	/*
	//test trace
	hmc_complex trace = trace_su3matrix(tmp2);
	trace.re -= 3.;
	complexaccumulate(&testsum,&trace);
	*/

	put_su3matrix(gaugefield,tmp2,spacepos,t,mu);
	/*	for(int n=0; n<2*SU3SIZE*VOLSPACE*NDIM*NTIME) {
	  gaugefield[n] = gf[n];
	  }*/
      }
    }
  }
  adjoin_su3(gaugefield,gaugefield2);
  adjoin_su3(gaugefield2,gaugefield);
  /*
  hmc_complex myctmp = global_trace_su3(gaugefield,2);
  testsum.re = myctmp.re - 3*VOLSPACE*NTIME;
  testsum.im = myctmp.im;
  */

  *check = testsum.re*testsum.re + testsum.im*testsum.im;
  *check += det_su3matrix(prod).re - hmc_one_f;
  }


  //CP: random number test
  int order[3]; 
  for(int i=0;i<size_1/3;i++){
  	random_1_2_3(order, &rnd[id]);
  	random_field_int[3*i] = order[0];
  	random_field_int[3*i + 1] = order[1];
  	random_field_int[3*i + 2] = order[2];
  }
  for(int i=0;i<size_2;i++){
    random_field_float[i] = ocl_new_ran( &rnd[id] );
  }
  hmc_float sample[4];
  
 
  SU2Update(sample, 3.555, &rnd[id]);
  
  su2mat[0] = sample[0];
  su2mat[1] = sample[1];
  su2mat[2] = sample[2];
  su2mat[3] = sample[3];
  
  //CP: heatbath test
  //take a link and its staplematrix after 200 host-iterations and run the host code and the device code on it
  hmc_ocl_su3matrix out_tmp[SU3SIZE];
  hmc_ocl_su3matrix in_tmp[SU3SIZE];
  hmc_ocl_staplematrix staple_tmp[STAPLEMATRIXSIZE];
  int heatbath_rnd_array_size = 10000;
  hmc_float heatbath_rnd_array[10000];
  for(int i = 0; i<heatbath_rnd_array_size; i++){
    heatbath_rnd_array[i] = heatbath_rnd_array_in[i];
  }
  
  for(int i = 0; i<SU3SIZE; i++){
    (in_tmp[i]).re = (heatbath_link_in[i]).re;
    (in_tmp[i]).im = (heatbath_link_in[i]).im;
  }
  for(int i = 0; i<STAPLEMATRIXSIZE; i++){
    (staple_tmp[i]).re = (heatbath_staple_in[i]).re;
    (staple_tmp[i]).im = (heatbath_staple_in[i]).im;
  }
  
  int cter;
  //CP: three different possibilities to test the heatbath
//   testing_heatbath_norandommat_no123(in_tmp, staple_tmp, out_tmp, 4.2);
//   testing_heatbath_no123(in_tmp, staple_tmp, out_tmp, 4.2, heatbath_rnd_array, &cter);
  testing_heatbath(in_tmp, staple_tmp, out_tmp, 4.2, heatbath_rnd_array, &cter);
  
  for(int i = 0; i<SU3SIZE; i++){
    (heatbath_link_out[i]).re = (out_tmp[i]).re;
    (heatbath_link_out[i]).im = (out_tmp[i]).im;
  }
  heatbath_cter[0] = cter;

  //CP: solver test: invert a small matrix on a cold-gaugeconfiguration and calculate the pion propagator

  simple_correlator(solver_spinor_in, solver_spinor_out, gaugefield2, solver_correlator, 0.125, 0.06, 0., 1000);

  return;
  } //else
}
