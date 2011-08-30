/** @file
 * Handling of lattice geometry.
 *
 * Utilities to map coordinates and indices as well as
 * move around in the index space.
 *
 * @todo As they are adresses, indices should use the type
 *       size_t instead of int.
 */
#ifndef _GEOMETRYH_
#define _GEOMETRYH_

#include "globaldefs.h"
#include <cstdlib>
#include <cmath>

//coord[0] = t
//coord[1] = x
//coord[2] = y
//coord[3] = z

/**
 * Calculate the index in even-odd preconditioned
 * storage based on the given non-preconditioned
 * indices.
 *
 * @param spacepos The index in the volume
 * @param timepos The index in time
 * @return even-odd preconditioned index
 */
int get_n_eoprec(int spacepos,int timepos);

//switch between (x,y,z) <-> nspace=0,...,VOLSPACE-1

/**
 * Calculate the spatial index of the given cartesian coordinates.
 *
 * @param coord Pointer to NDIM integers representing cartesian coordinates.
 *              time is expected in index 0 of this array.
 * @param Spatial index
 */
int get_nspace(int* coord);
/**
 * Calculate the cartesian coordinates of the given spatial index.
 *
 * @param nspace Spatial index
 * @param dir The dimension for which to retrieve the cartisian coordinate
 * @param Cartesian coordinate of nspace in dimension dir.
 */
int get_spacecoord(int nspace, int dir);

//get spatial neighbors

/**
 * Calculate the next spatial index in a given direction.
 *
 * @param nspace The spatial index to start with
 * @param dir The dimension in which to move
 * @return Index of the next site in the given direction
 */
int get_neighbor(int nspace, int dir);
/**
 * Calculate the previous spatial index in a given direction.
 *
 * @param nspace The spatial index to start with
 * @param dir The dimension in which to move
 * @return Index of the next site in the given direction
 */
int get_lower_neighbor(int nspace, int dir);
//Checkerboard: get real coordinates from EVEN/ODD-index
/**
 * Get the non-preconditioned indices for a given even-odd preconditioned index
 * of an even site.
 *
 * @param[in] idx The index
 * @param[out] out_space Pointer to which to write the spatial index.
 * @param[out] out_t Pointer to which to write the index in time direction.
 */
void get_even_site(int idx, int * out_space, int * out_t);
/**
 * Get the non-preconditioned indices for a given even-odd preconditioned index
 * of an odd site.
 *
 * @param[in] idx The index
 * @param[out] out_space Pointer to which to write the spatial index.
 * @param[out] out_t Pointer to which to write the index in time direction.
 */
void get_odd_site(int idx, int * out_space, int * out_t);

/**
 * Get the non-even-odd-preconditioned index of a site based on the spatial and temporal
 * index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_pos(int spacepos, int t);

/**
 * Get the non-even-odd-preconditioned index link based on the spatial, temporal
 * and Dirac index.
 *
 * @param spacepos Spatial index
 * @param t Temporal index
 * @return Global index
 */
int get_global_link_pos(int mu, int spacepos, int t);

//get gaugefield element from long array
/** This function returns the index of a specific hmc_float out of an hmc_complex array. 
 *  This is used to when copying the gaugefield to the OpenCL device.
 *  One has:
 *   c: complex index (0 for real, 1 for imaginary part)
 *   a, b: SU3 indices
 *   mu, spacepos, t: Dirac and spacetime indices
 * @todo this function is highly misterious. Still??
 */
int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t);

//Spinor functions

/**
 * Get the index of a spinor element within the spinor
 *
 * @param alpha Spin index
 * @param color Color index
 * @return Spinor local index
 */
int spinor_element(int alpha, int color);
/**
 * Get the index of a spinor element within the non-even-odd-preconditioned lattice.
 *
 * @param alpha Spin index
 * @param color Color index
 * @param nspace Spatial index
 * @param t Temporal index
 * @return Element index within the lattice
 */
int spinor_field_element(int alpha, int color, int nspace, int t);
/**
 * Get the color component from a spinor-local index.
 *
 * @param spinor_element Spinor-local indexx
 * @return Color component
 */
int spinor_color(int spinor_element);
/**
 * Get the spint component from a spinor-local index.
 *
 * @param spinor_element Spinor-local indexx
 * @return Spin component
 */
int spinor_spin(int spinor_element,int color);
/**
 * Get the index of a spinor element within the even-odd-preconditioned lattice.
 *
 * @param alpha Spin index
 * @param color Color index
 * @param nspace Spatial index
 * @param t Temporal index
 * @return Element index within the lattice
 */
int eoprec_spinor_field_element(int alpha, int color, int nspace, int t);
/**
 * Get the index of a spinor element within the even-odd-preconditioned lattice.
 *
 * @param alpha Spin index
 * @param color Color index
 * @param n_eoprec Even-odd-preconditioned index as returned by get_n_eoprec(int,int).
 * @return Element index within the lattice
 */
int eoprec_spinor_field_element(int alpha, int color, int n_eoprec);

#endif
