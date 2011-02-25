#ifndef _GAUGEOBSERVABLESH_
#define _GAUGEOBSERVABLESH_
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_geometry.h"
#include "host_use_timer.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <cmath>

void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2);
void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2, int iter);
void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2, int iter, std::string file);
void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter);
void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter, std::string file);

hmc_float plaquette(hmc_gaugefield * field, hmc_float* tplaq, hmc_float* splaq);
hmc_float plaquette(hmc_gaugefield * field);
hmc_complex polyakov(hmc_gaugefield * field);
hmc_complex spatial_polyakov(hmc_gaugefield * field, int dir);

#endif
