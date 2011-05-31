#include "host_solver.h"

hmc_error solver_eoprec(inputparameters * parameters, hmc_spinor_field* in,  hmc_eoprec_spinor_field* be, hmc_eoprec_spinor_field* bo, hmc_gaugefield* gaugefield, int use_cg, hmc_spinor_field* out){	
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
	
	
  //!!CP: in the end this should be done with a compiler option
  use_cg = FALSE;
  hmc_eoprec_spinor_field even[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field odd[EOPREC_SPINORFIELDSIZE];
	hmc_eoprec_spinor_field tmp[EOPREC_SPINORFIELDSIZE];

	convert_to_eoprec(even,odd,in);
  convert_to_kappa_format_eoprec(even,( (*parameters).get_kappa() ));
  convert_to_kappa_format_eoprec(odd,( (*parameters).get_kappa() ));
		
	//calculate even-solution even
	if(use_cg)
		cg_eoprec(parameters, even, be, gaugefield);
	else
		bicgstab_eoprec(parameters, even, be, gaugefield);
		
	//desired solution is g = L^-1 x
	//calculate L^-1x_even = g_even
	hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
	dslash_eoprec(parameters, even,gaugefield, ODD , tmp); //use newood as tmp variable
	M_inverse_sitediagonal(parameters, tmp, spintmp);
		
	//calculate odd-solution odd
	M_inverse_sitediagonal(parameters, bo, odd);
		
	//calculate g_odd
	hmc_complex one = hmc_complex_one;
	saxpy_eoprec(spintmp, odd, &one, odd);

	convert_from_kappa_format_eoprec(even,even, ( (*parameters).get_kappa() ));
	convert_from_kappa_format_eoprec(odd, odd, ( (*parameters).get_kappa() ));

	//write out g
	convert_from_eoprec(even,odd,out);
  
	return HMC_SUCCESS;
}

//here, one can switch with use_cg = TRUE/FALSE which Matrix should be inverted
//TODO think about if the use of CG = MdaggerM and BiCGStab = M is the best solution
hmc_error solver(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* b, hmc_gaugefield* gaugefield, int use_cg, hmc_spinor_field* out){
	
	convert_to_kappa_format(in, ( (*parameters).get_kappa() ));
	int err = 0;
	if(use_cg)
		err = cg(parameters, in, b, gaugefield);
	else
		err = bicgstab(parameters, in, b, gaugefield);
	convert_from_kappa_format(in, out, ( (*parameters).get_kappa() ));
	if(err != HMC_SUCCESS) return HMC_STDERR;
	    
	return HMC_SUCCESS;
}

hmc_error bicgstab(inputparameters * parameters, hmc_spinor_field* inout, hmc_spinor_field* source, hmc_gaugefield* gaugefield){

  //BiCGStab according to hep-lat/9404013

	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();


  hmc_spinor_field* rn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* rhat = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* v = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* p = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* s = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* t = new hmc_spinor_field[SPINORFIELDSIZE];
	hmc_spinor_field* aux = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;

  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      M(parameters, inout,gaugefield,rn);
      saxpy(rn, source, &one, rn);
      copy_spinor(rn, rhat);

      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;
      set_zero_spinorfield(v);
      set_zero_spinorfield(p);
      
      //printf("true residue squared: %e\n",global_squarenorm(rn));
    }

    rho_next = scalar_product(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz(p, v, rn, &beta, &tmp2, p);

    M(parameters, p,gaugefield,v);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy(v, rn, &alpha, s);

    M(parameters, s,gaugefield,t);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);
    omega = complexdivide(&(tmp1),&(tmp2));

    saxpy(t, s, &omega, rn);

    saxsbypz(p, s, inout, &alpha, &omega, inout);

    hmc_float resid = global_squarenorm(rn);

    if(resid<epssquare) {
      M(parameters, inout,gaugefield,aux);
      saxpy(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm(aux);
      //printf("%d\t%e\t%e\n",iter,resid,trueresid);
      if(trueresid<epssquare) return HMC_SUCCESS;
    }
  }

  delete [] rn;
  delete [] rhat;
  delete [] v;
  delete [] p;
  delete [] s;
  delete [] t;
	delete [] aux;
	
  return HMC_SUCCESS;
}


hmc_error bicgstab_eoprec(inputparameters * parameters, hmc_eoprec_spinor_field* inout,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield){

	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
	
  //BiCGStab according to hep-lat/9404013
  hmc_eoprec_spinor_field* rn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* rhat = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* v = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* p = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* s = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* t = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* aux = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];

	hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;

  for(int iter=0; iter<cgmax; iter++){
    if(iter%iter_refresh==0) {
      //fresh start
			set_zero_spinorfield_eoprec(p);
      Aee(parameters, inout,gaugefield,rn);
      saxpy_eoprec(rn, source, &one, rn);
      copy_spinor_eoprec(rn, rhat);
      
      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;      
      set_zero_spinorfield_eoprec(v);
      
      
//       printf("true residue squared: %e\n",global_squarenorm_eoprec(rn));
    }

    rho_next = scalar_product_eoprec(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz_eoprec(p, v, rn, &beta, &tmp2, p);
    
    Aee(parameters, p,gaugefield,v);

    tmp1 = scalar_product_eoprec(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy_eoprec(v, rn, &alpha, s);

    Aee(parameters, s,gaugefield,t);

    tmp1 = scalar_product_eoprec(t,s);
    tmp2 = scalar_product_eoprec(t,t);
    omega = complexdivide(&tmp1,&tmp2);

    saxpy_eoprec(t, s, &omega, rn);

    saxsbypz_eoprec(p, s, inout, &alpha, &omega, inout);
  
    hmc_float resid = global_squarenorm_eoprec(rn);
//     printf("%d\t%e\n",iter,resid);

    if(resid<epssquare) {
      Aee(parameters, inout,gaugefield,aux);
      saxpy_eoprec(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm_eoprec(aux);
//       printf("true residue squared: %e\n",trueresid);
      if(trueresid<epssquare) {
				return HMC_SUCCESS;     
			}
    }
  }

  delete [] rn;
  delete [] rhat;
  delete [] v;
  delete [] p;
  delete [] s;
  delete [] t;
	delete [] aux;
	
  return HMC_SUCCESS;
}

//CP: this is defined directly with QplusQminus
hmc_error cg(inputparameters * parameters, hmc_spinor_field* inout, hmc_spinor_field* source, hmc_gaugefield* gaugefield){
 
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
	
	
	hmc_spinor_field* rn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* xnn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* pn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* pnn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* rnn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* tmp = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_complex alpha;
	hmc_complex beta;
	hmc_complex tmp1;
	hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;
	
	int iter;
	//debugging
	int cter =0;
  //main loop
  for(iter = 0; iter < cgmax; iter ++){  
    if(iter%iter_refresh==0){
			//fresh start
			QplusQminus(parameters, inout,gaugefield,rn);
			saxpy(rn, source, &one, rn);
			copy_spinor(rn, pn);
// 			printf("true residue squared: %e\n",global_squarenorm(rn));
		}
		
		QplusQminus(parameters, pn,gaugefield,tmp);
		tmp1 = scalar_product(rn, pn);
		tmp2 = scalar_product(pn, tmp);
		alpha = complexdivide(&tmp1, &tmp2);
		
		saxpy(tmp, rn, &alpha, rnn);
		alpha = complexmult(&alpha, &minusone);
		saxpy(pn, inout, &alpha, xnn);
		
		hmc_float resid = global_squarenorm(rnn);
    
    if (resid < epssquare){
// 			printf("residue small enough: %e\n",resid);
			copy_spinor(xnn, inout);
			break;
    }
    else{
// 			printf("residue not small enough: %e\n",resid);
			tmp1 = scalar_product(rn, rn);
			tmp2 = scalar_product(rnn, rnn);
			beta = complexdivide(&tmp2, &tmp1);
			beta = complexmult(&beta, &minusone);
			saxpy(pn, rnn, &beta, pnn);
			//copy spinors for next iteration
			copy_spinor(pnn, pn);
			copy_spinor(rnn, rn);
			copy_spinor(xnn, inout);
			
			//debugging
			cter++;
    }
  }
	
	delete [] rn;
	delete [] xnn;
	delete [] pn;
	delete [] pnn;
	delete [] rnn;
	delete [] tmp;
// 	printf("cter: %i \n\n", cter);
	if(cter < cgmax) return HMC_SUCCESS;
	if(cter == cgmax) return HMC_STDERR;
}

/** @todo CP: this cannot be used with Aee since it is not hermitian, instead one has to insert the eoprec-version of QplusQminus!!! */
hmc_error cg_eoprec(inputparameters * parameters, hmc_eoprec_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield){
 
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
	
	
	hmc_eoprec_spinor_field* rn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* xnn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* pn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* pnn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* rnn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* tmp = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_complex alpha;
	hmc_complex beta;
	hmc_complex tmp1;
	hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;
	
	int iter;

  for(iter = 0; iter < cgmax; iter ++){  
    if(iter%iter_refresh==0){
			//fresh start
			Aee(parameters, inout,gaugefield,rn);
			saxpy_eoprec(rn, source, &one, rn);
			copy_spinor_eoprec(rn, pn);
// 			printf("true residue squared: %e\n",global_squarenorm(rn));
		}
		
		Aee(parameters, pn,gaugefield,tmp);
		tmp1 = scalar_product_eoprec(rn, pn);
		tmp2 = scalar_product_eoprec(pn, tmp);
		alpha = complexdivide(&tmp1, &tmp2);
		
		saxpy_eoprec(tmp, rn, &alpha, rnn);
		alpha = complexmult(&alpha, &minusone);
		saxpy_eoprec(pn, inout, &alpha, xnn);
		
		hmc_float resid = global_squarenorm_eoprec(rnn);
    
    if (resid < epssquare){
			printf("residue small enough: %e\n",resid);
			copy_spinor_eoprec(xnn, inout);
			break;
    }
    else{
			tmp1 = scalar_product_eoprec(rn, rn);
			tmp2 = scalar_product_eoprec(rnn, rnn);
			beta = complexdivide(&tmp2, &tmp1);
			beta = complexmult(&beta, &minusone);
			saxpy_eoprec(pn, rnn, &beta, pnn);
			//copy spinors for next iteration
			copy_spinor_eoprec(pnn, pn);
			copy_spinor_eoprec(rnn, rn);
			copy_spinor_eoprec(xnn, inout);
    }
  }
	
	delete [] rn;
	delete [] xnn;
	delete [] pn;
	delete [] pnn;
	delete [] rnn;
	delete [] tmp;
	
	return HMC_SUCCESS;
}

