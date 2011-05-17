/** @file
 * Observables calculated from the gauge field
 */

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

/**
 * Print calculated plaquette and polyakov to stdout.
 *
 * @param field The gaugefield from which to calculate the observables
 * @param timer The timer in which to aggregate plaquette calculation time
 * @param timer2 The timer in which to aggregate plyakov calculation time
 */
void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2);
/**
 * Print iteration and calculated plaquette and polyakov to stdout.
 *
 * @param field The gaugefield from which to calculate the observables
 * @param timer The timer in which to aggregate plaquette calculation time
 * @param timer2 The timer in which to aggregate plyakov calculation time
 * @param iter The number to use for the iteration
 */
void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2, int iter);
/**
 * Print iteration and calcualted plaquette and polyakov to a file.
 *
 * @param field The gaugefield from which to calculate the observables
 * @param timer The timer in which to aggregate plaquette calculation time
 * @param timer2 The timer in which to aggregate plyakov calculation time
 * @param iter The number to use for the iteration
 * @param file The name of the file to print to
 */
void print_gaugeobservables(hmc_gaugefield* field, usetimer* timer, usetimer * timer2, int iter, std::string file);
/**
 * Print the passed gaugeobservable to stdout.
 */
void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter);
/**
 * Print the passed gaugeobservable to a file.
 */
void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter, std::string file);

/**
 * Calculate the plaquette from the given gaugefield.
 */
hmc_float plaquette(hmc_gaugefield * field, hmc_float* tplaq, hmc_float* splaq);
/**
 * Calculate the plaquette from the given gaugefield.
 */
hmc_float plaquette(hmc_gaugefield * field);
/**
 * Calculate the polyakov from the given gaugefield.
 */
hmc_complex polyakov(hmc_gaugefield * field);
/**
 * Calculate the spatial polyakov from the given gaugefield.
 */
hmc_complex spatial_polyakov(hmc_gaugefield * field, int dir);

#endif /* _GAUGEOBSERVABLESH_ */
