#include "host_gaugefieldoperations.h"

hmc_error init_gaugefield(hmc_gaugefield* gaugefield,inputparameters* parameters, usetimer* timer){
  sourcefileparameters parameters_source;
  if((*parameters).get_startcondition()==START_FROM_SOURCE){
    int err;
    (*timer).reset();
    //tmp gauge field
    hmc_float * gaugefield_tmp;
    gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float)*NDIM*NC*NC*NTIME*VOLSPACE);
    err = parameters_source.readsourcefile(&((*parameters).sourcefile)[0], (*parameters).get_prec(), &gaugefield_tmp);
    err = copy_gaugefield_from_ildg_format(gaugefield, gaugefield_tmp, parameters_source.num_entries_source);
    free(gaugefield_tmp);
    (*timer).add();
    if (err == 0){
      print_info_source(&parameters_source);
    } else {
      printf("error in setting vals from source!!! check global settings!!!\n\n");
      return HMC_XMLERROR;
    }
  }
  if((*parameters).get_startcondition()==COLD_START) {
    (*timer).reset();
    set_gaugefield_cold(gaugefield);
    (*timer).add();
  }
  if((*parameters).get_startcondition()==HOT_START) {
    (*timer).reset();
    set_gaugefield_cold(gaugefield);
    (*timer).add();
  }
  return HMC_SUCCESS;
}

void print_info_source(sourcefileparameters* params){
  printf("**********************************************************\n");
  printf("Sourcefile parameters: (list not complete)\n");
  printf("field:  %s\n",(*params).field_source);
  printf("LX:  \t%d\n",(*params).lx_source);
  printf("LY:  \t%d\n",(*params).ly_source);
  printf("LZ:  \t%d\n",(*params).lz_source);
  printf("LT:  \t%d\n",(*params).lt_source);
  printf("entries: %d\n", (*params).num_entries_source);
  printf("beta:  \t%f\n",(*params).beta_source);
  printf("mu:  \t%f\n",(*params).mu_source);
  printf("kappa:  %f\n",(*params).kappa_source);
  printf("mubar:  %f\n",(*params).mubar_source);
  printf("plaq: \t%f\n",(*params).plaquettevalue_source);
  printf("**********************************************************\n");
  printf("\n");
  return;
}
