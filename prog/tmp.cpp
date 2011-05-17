/* Updates the momenta: equation 16 of Gottlieb */
void update_momenta(int * mnllist, double step, const int no) {
  int i,mu, k;
  double tmp;
  su3adj *xm,*deriv;
  double sum=0., max=0.;
  double sum2=0.;
  double atime=0., etime=0.;

  for(i=0;i<(VOLUME);i++) { 
    for(mu=0;mu<4;mu++) { 
//df0 is the global momentum array
      _zero_su3adj(df0[i][mu]);
    }
  }

  for(k = 0; k < no; k++) {
    sum = 0.;
    max = 0.;
    for(i = (VOLUME); i < (VOLUME+RAND); i++) { 
      for(mu = 0; mu < 4; mu++) { 
	_zero_su3adj(df0[i][mu]);
      }
    }
    /* the next line is the only part which will change between PBC and SF:
     when the "GAUGE MONOMIAL / SFGAUGE MONOMIAL" is chosen
     the correct derivative function for each case is used "gauge_derivative / sf_gauge_derivative" */
    monomial_list[ mnllist[k] ].derivativefunction(mnllist[k]);
    for(i = 0; i < VOLUME; i++) {
      for(mu = 0; mu < 4; mu++) {
	xm=&moment[i][mu];
	deriv=&df0[i][mu];
	/* force monitoring */
	if(g_debug_level > 0) {
	  sum2 = _su3adj_square_norm(*deriv); 
	  sum+= sum2;
	  if(fabs(sum2) > max) max = sum2;
	}
	tmp = step*monomial_list[ mnllist[k] ].forcefactor;
	/* the minus comes from an extra minus in trace_lambda */
	_minus_const_times_mom(*xm,tmp,*deriv); 
	/* set to zero immediately */
	_zero_su3adj(df0[i][mu]);
      }
    }
  }
  return;
}