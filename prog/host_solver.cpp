#include "host_solver.h"

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){
  //CP: in the end this should be done with a compiler option
  int use_cg = 0;
  if(!use_eo){
    hmc_eoprec_spinor_field even[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field odd[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field newodd[EOPREC_SPINORFIELDSIZE];
		hmc_eoprec_spinor_field tmp[EOPREC_SPINORFIELDSIZE];

    hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
    
		convert_to_eoprec(even,odd,in);
    convert_to_kappa_format_eoprec(even,kappa);
    convert_to_kappa_format_eoprec(odd,kappa);
 
    create_point_source_eoprec(be,bo,1,0,0,kappa,mu,theta, chem_pot_re, chem_pot_im, gaugefield);
		
		//calculate even-solution even
		if(use_cg)
			cg_eoprec(even, be, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);
		else
    bicgstab_eoprec(even, be, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);
		
		//desired solution is g = L^-1 x
		//calculate L^-1x_even = g_even
    hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
    dslash_eoprec(even, tmp,gaugefield,kappa, theta, chem_pot_re, chem_pot_im, ODD); //use newood as tmp variable
    M_inverse_sitediagonal(tmp, spintmp,kappa,mu);
		
		//calculate odd-solution odd
    M_inverse_sitediagonal(bo, odd, kappa,mu);
		
		//calculate g_odd
		hmc_complex one = hmc_complex_one;
		saxpy_eoprec(spintmp, odd, &one, odd);

    convert_from_kappa_format_eoprec(even,even, kappa);
    convert_from_kappa_format_eoprec(odd, odd, kappa);
		//write out g
    convert_from_eoprec(even,odd,out);
  
    return HMC_SUCCESS;
  }
  else{
    convert_to_kappa_format(in, kappa);
		if(use_cg)
			cg(in, b, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);
		else
			bicgstab(in, b, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);
    convert_from_kappa_format(in, out, kappa);
    
    return HMC_SUCCESS;
  }
}

hmc_error bicgstab(hmc_spinor_field* inout, hmc_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){

  //BiCGStab according to hep-lat/9404013

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
      M(inout,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
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

    M(p,v,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy(v, rn, &alpha, s);

    M(s,t,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);
    omega = complexdivide(&(tmp1),&(tmp2));

    saxpy(t, s, &omega, rn);

    saxsbypz(p, s, inout, &alpha, &omega, inout);

    hmc_float resid = global_squarenorm(rn);

    if(resid<epssquare) {
      M(inout,aux,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
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


hmc_error bicgstab_eoprec(hmc_eoprec_spinor_field* inout,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im,int cgmax){

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
      Aee(inout,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy_eoprec(rn, source, &one, rn);
      copy_spinor_eoprec(rn, rhat);
      
      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;      
      set_zero_spinorfield_eoprec(v);
      set_zero_spinorfield_eoprec(p);
      
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
    
    Aee(p,v,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product_eoprec(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy_eoprec(v, rn, &alpha, s);

    Aee(s,t,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product_eoprec(t,s);
    tmp2 = scalar_product_eoprec(t,t);
    omega = complexdivide(&tmp1,&tmp2);

    saxpy_eoprec(t, s, &omega, rn);

    saxsbypz_eoprec(p, s, inout, &alpha, &omega, inout);
  
    hmc_float resid = global_squarenorm_eoprec(rn);
//     printf("%d\t%e\n",iter,resid);

    if(resid<epssquare) {
      Aee(inout,aux,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy_eoprec(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm_eoprec(aux);
//       printf("true residue squared: %e\n",trueresid);
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

hmc_error cg(hmc_spinor_field* inout, hmc_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){
 
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
  //main loop

  for(iter = 0; iter < cgmax; iter ++){  
    if(iter%iter_refresh==0){
			//fresh start
			M(inout,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
			saxpy(rn, source, &one, rn);
			copy_spinor(rn, pn);
			printf("true residue squared: %e\n",global_squarenorm(rn));
		}
		
		M(pn,tmp,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
		tmp1 = scalar_product(rn, pn);
		tmp2 = scalar_product(pn, tmp);
		alpha = complexdivide(&tmp1, &tmp2);
		
		saxpy(tmp, rn, &alpha, rnn);
		alpha = complexmult(&alpha, &minusone);
		saxpy(pn, inout, &alpha, xnn);
		
		hmc_float resid = global_squarenorm(rnn);
    
    if (resid < epssquare){
			printf("residue small enough: %e\n",resid);
			copy_spinor(xnn, inout);
			break;
    }
    else{
			tmp1 = scalar_product(rn, rn);
			tmp2 = scalar_product(rnn, rnn);
			beta = complexdivide(&tmp2, &tmp1);
			beta = complexmult(&beta, &minusone);
			saxpy(pn, rnn, &beta, pnn);
			//copy spinors for next iteration
			copy_spinor(pnn, pn);
			copy_spinor(rnn, rn);
			copy_spinor(xnn, inout);
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

hmc_error cg_eoprec(hmc_eoprec_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){
 
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
			Aee(inout,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
			saxpy_eoprec(rn, source, &one, rn);
			copy_spinor_eoprec(rn, pn);
// 			printf("true residue squared: %e\n",global_squarenorm(rn));
		}
		
		Aee(pn,tmp,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
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

