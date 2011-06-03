/** @file
 * Mathematical operations on the fermion matrix
 */

#ifndef _OPERATIONS_FERMIONMATRIXH_
#define _OPERATIONS_FERMIONMATRIXH_

#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_solver.h"
#include <cmath>

/**
 * Calculates $f/ Q^+Q^⁻ \psi = \left( \gamma_5 M_{\tilde{\mu}}} \gamma_5 M_{-\tilde{\mu}}}\right ) \psi /f$
 * @param[in] parameters includes parameters needed for evaluation of M
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 */
hmc_error QplusQminus(inputparameters * parameters, hmc_spinor_field * in, hmc_gaugefield * field, hmc_spinor_field * out);

/**
 * Calculates $f/ Q^+ \psi = \left( \gamma_5 M_{\tilde{\mu}}} \right ) \psi /f$
 * @param[in] parameters includes parameters needed for evaluation of M
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 */
hmc_error Qplus(inputparameters * parameters, hmc_spinor_field * in, hmc_gaugefield * field, hmc_spinor_field * out);

/**
 * Calculates $f/ Q^- \psi = \left( \gamma_5 M_{-\tilde{\mu}}} \right ) \psi /f$
 * @param[in] parameters includes parameters needed for evaluation of M
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 */
hmc_error Qminus(inputparameters * parameters, hmc_spinor_field * in, hmc_gaugefield * field, hmc_spinor_field * out);

/**
 * Calculates $f/ M \psi = \left( M_{diag} - \kappa\notD\right ) \psi /f$
 * @param[in] parameters includes parameters needed for evaluation of M
 * @param[in] in spinorfield that M acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error M(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out);

/**
 * Calculates $f/ \notD \psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error dslash(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out);

/**
 * Calculates $f/M_{diag}  \psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error M_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out);//CP: checked

/**
 * Calculates $f/ \notD \psi /f$ in the eoprec-version
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[in] evenodd act on even or odd site
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error dslash_eoprec(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, int evenodd, hmc_eoprec_spinor_field* out);

/**
 * Calculates $f/M_{diag}^{-1}  \psi /f$ for eoprec fields
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error M_inverse_sitediagonal(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out);//CP: checked

/**
 * Calculates $f/M_{diag}  \psi /f$ for eoprec fields
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error M_sitediagonal(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out);//CP: checked

/**
 * Calculates $A_{ee} \psi /f$ for eoprec fields
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 */
hmc_error Aee(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_eoprec_spinor_field* out);

//helper functions

/**
 * Calculates $/f \not D \psi(specific site) = ... /f$ in temporal direction ($f/ \mu = 0 /f$)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites 
 * @param[out] spinout output spinor
 * @remark tested by CP
 */
void dslash_temporal (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked

/**
 * Calculates $/f \not D \psi(specific site) = ... /f$ in spatial direction ($f/ \mu = 1,2,3 /f$)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark tested by CP
 */
void dslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked

/**
 * Calculates $/f \not D \psi(specific site) = ... /f$ in temporal direction ($f/ \mu = 0 /f$) (EOPREC-version)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark tested by CP
 */
void dslash_temporal_eoprec (hmc_spinor * spinout, int pos, int t, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked

/**
 * Calculates $/f \not D \psi(specific site) = ... /f$ in spatial direction ($f/ \mu = 1,2,3 /f$) (EOPREC-version)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark tested by CP
 */
void dslash_spatial_eoprec (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecated functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//CP: this was usefull because M_daggerM can also be calculated by successivly using M and gamma_5

/**
 * Calculates $f/ M^{\dagger} \psi = \left( M_{diag}^{\dagger} + \notD^{\dagger} \right ) \psi /f$
 * @param[in] parameters includes parameters needed for evaluation of $f/ M^{\dagger} /f$
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 * @remark deprecated by CP
 */
// hmc_error Mdagger(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out);

/**
 * Calculates $f/ \notD\dagger} \psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 * @remark deprecated by CP
 */
// hmc_error ddaggerslash(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out);

/**
 * Calculates $f/M_{diag}^{\dagger}  \psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[out] out output spinorfield
 * @remark deprecated by CP
 */
// hmc_error Mdagger_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out);//CP: not checked

/**
 * Calculates $f/ M^{\dagger} \psi = \left( M_{diag} + \notD\right ) \left( M_{diag} + \notD\right )^{\dagger}\psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark tested by CP
 * @remark deprecated by CP
 */
// hmc_error MdaggerM(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out);

/**
 * Calculates $f/ \left( M_{diag}\right ) \left( M_{diag}  )^{\dagger}\psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[out] out output spinorfield
 * @remark deprecated by CP
 */
// hmc_error MdaggerM_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out); //CP: not checked

/**
 * Calculates $f/  \left( \notD\right ) \left( \notD\right )^{\dagger}\psi /f$
 * @param[in] parameters includes parameters needed
 * @param[in] in spinorfield that the matrix acts on
 * @param[in] gaugefield input gaugefield
 * @param[out] out output spinorfield
 * @remark deprecated by CP
 */
// hmc_error ddaggerd(inputparameters * parameter, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out); //CP: not checked


/**
 * Calculates $/f \not D^{\dagger} \psi(specific site) = ... /f$ in temporal direction ($f/ \mu = 0 /f$)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark tested by CP
 * @remark deprecated by CP
 */
// void ddaggerslash_temporal (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked

/**
 * Calculates $/f \not D^{\dagger} \psi(specific site) = ... /f$ in spatial direction ($f/ \mu = 1,2,3 /f$)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark tested by CP
 * @remark deprecated by CP
 */
// void ddaggerslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);//CP: checked



/**
 * Calculates $/f \not D^{\dagger} \not D \psi(specific site) = ... /f$ in all directions ($f/ \mu = 0,1,2,3 /f$)
 * @param[in] parameters includes parameters needed
 * @param[in] pos spatial position of site
 * @param[in] t temporal position of site
 * @param[in] in inputspinorfield that include the spinors acting on the one at site (pos, t)
 * @param[in] gaugefield input gaugefield providing the links between sites
 * @param[out] spinout output spinor
 * @remark deprecated by CP
 */
// void ddaggerd_calc (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im);



#endif /* _OPERATIONS_FERMIONMATRIXH_ */
