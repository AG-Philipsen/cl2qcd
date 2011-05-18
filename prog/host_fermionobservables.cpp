#include "host_fermionobservables.h"

hmc_error simple_correlator(inputparameters * parameters, hmc_gaugefield* gaugefield){

	//CP: one needs bicgstab here for M
	int use_cg = FALSE;
	
  hmc_spinor_field in[SPINORFIELDSIZE];
	hmc_spinor_field phi[SPINORFIELDSIZE];
	init_spinorfield_cold(in);

  //pseudo scalar, flavour multiplet
  hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }

  for(int k=0; k<NC*NSPIN; k++) {
		if(!use_eo){
			hmc_spinor_field b[SPINORFIELDSIZE];
			create_point_source(parameters,k,0,0,b);
			solver(parameters, in, b, gaugefield, use_cg, phi);
		}
		else{
			hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
			hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
			
			create_point_source_eoprec(parameters, k,0,0, gaugefield, be,bo);
			solver_eoprec(parameters, in, be, bo, gaugefield, use_cg, phi);
			
		}
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

  printf("pseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }
  return HMC_SUCCESS;
}
