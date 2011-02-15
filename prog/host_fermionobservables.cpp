#include "host_fermionobservables.h"

hmc_error simple_correlator(hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  hmc_spinor_field in[SPINORFIELDSIZE];
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(in, n, t, a, j);
	}
      }
    }
  }
  hmc_float norm=global_squarenorm(in);
  norm = sqrt(norm);
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    in[n].re /= norm;
    in[n].im /=norm;
  }


  //pseudo scalar, flavour multiplet
  hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }

  hmc_spinor_field b[SPINORFIELDSIZE];
  hmc_spinor_field phi[SPINORFIELDSIZE];

  for(int k=0; k<NC*NSPIN; k++) {
    create_point_source(b,k,0,0,kappa,mu,gaugefield);
    solver(in, phi, b, gaugefield, kappa, mu, theta, cgmax);

    for(int timepos = 0; timepos<NTIME; timepos++) {
      for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
	for(int alpha = 0; alpha<NSPIN; alpha++) {
	  for(int c = 0; c<NC; c++) {
	    //	    int j = spinor_element(alpha,c);
	    int n = spinor_field_element(alpha, c, spacepos, timepos);
	    int z = get_spacecoord(spacepos, 3);
	    hmc_complex tmp = phi[n];
	    hmc_complex ctmp = complexconj(&tmp);
	    hmc_complex incr = complexmult(&ctmp,&tmp);
	    correlator_ps[z].re += incr.re;
	    correlator_ps[z].im += incr.im;
	  }
	}
      }
    }
  }

  printf("pseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }
  return HMC_SUCCESS;
}
