//opencl_geometry.cl
/** @file
 * Device code for lattice geometry handling
 * Geometric conventions used should only be applied in this file
 * All external functions should only use functions herein!!!
 */

////////////////////////////////////////////////////////////////
// General Definitions
////////////////////////////////////////////////////////////////

/** In the following, identify (coord being coord_full or coord_spatial):
 *coord.x == x 
 *coord.y == y 
 *coord.z == z 
 *(coord.w == t)
 * NOTE: This does not necessarily reflect the geometric conventions used! 
 */

/** index type to store all spacetime coordinates  */
typedef uint4 coord_full;

/** index type to store all spatial coordinates */
typedef uint3 coord_spatial;

/** index type to store a temporal coordinate */ 
typedef uint coord_temporal;

/** index type to store a site/link position */ 
typedef uint site_idx;
typedef uint link_idx;

/** index type to store a spatial site position */ 
typedef uint spatial_idx;

/** index type to store a direction (0..3) */
typedef uint dir_idx;

/** An index type storing space (index) and time (coordinate) seperate. */
typedef struct {
  spatial_idx space;
  coord_temporal time;
} st_idx;

////////////////////////////////////////////////////////////////
// Geometric Conventions
////////////////////////////////////////////////////////////////

/** Identify each spacetime direction */
#define TDIR 0
#define XDIR 1
#define YDIR 2
#define ZDIR 3

/**
 * The following conventions are used:
 * (NS: spatial extent, NT: temporal extent, NDIM: # directions, VOL4D: lattice volume, VOLSPACE: spatial volume)
 * A spatial idx is adressed as spatial_idx(x,y,z) = x * NS^(XDIR-1) + y * NS^(YDIR-1) + z * NS^(ZDIR-1)
 * A site idx is addressed as site_idx(x,y,z,t) = spatial_idx(x,y,z) + t*NS*NS*NS
 * A link idx is addressed as link_idx(x,y,z,t,mu) = mu + NDIM*site_idx(x,y,z,t)
 */

/**
 * The following conventions are used in tmlqcd:
 * #define TDIR 0
 * #define XDIR 3
 * #define YDIR 2
 * #define ZDIR 1
 * (NS: spatial extent, NT: temporal extent, NDIM: # directions, VOL4D: lattice volume, VOLSPACE: spatial volume)
 * A spatial idx is adressed as spatial_idx(x,y,z) = x * NS^(XDIR-1) + y * NS^(YDIR-1) + z * NS^(ZDIR-1) (same)
 * A site idx is addressed as site_idx(x,y,z,t) = spatial_idx(x,y,z) + t * VOLSPACE (same)
 * A link idx is addressed as an array with 2 indices.
 */

////////////////////////////////////////////////////////////////
// Functions relying explicitely on the geometric conventions
// defined above (w/o even-odd functions!!)
////////////////////////////////////////////////////////////////

/** spatial coordinates <-> spatial_idx
 *using the convention:
 *spatial_idx = x + NS * y + NS*NS * z
 *@todo this can be generalize using the definitions of the spatial directions
 *see  http://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int 
 */
// spatial_idx get_nspace(const coord_spatial coord)
spatial_idx get_spatial_idx(const coord_spatial coord)
{
  return (coord.x +  NSPACE * coord.y + NSPACE * NSPACE * coord.z);
  // return (coord.x * NSPACE^(XDIR-1) + coord.y * NSPACE^(YDIR-1) + coord.z * NSPACE^(ZDIR-1) );
}
coord_spatial get_coord_spatial(const spatial_idx nspace)
{
  coord_spatial coord;
  coord.z = nspace / NSPACE / NSPACE;
  uint acc = coord.z;
  coord.y = nspace / NSPACE - NSPACE * acc;
  acc = NSPACE * acc + coord.y;
  coord.x = nspace - NSPACE * acc;
  return coord;
}

/** 
 * st_idx <-> site_idx using the convention:
 * site_idx = x + y*NS + z*NS*NS + t*NS*NS*NS
 * = site_idx_spatial + t*VOLSPACE
 * <=>t = site_idx / VOLSPACE
 * site_idx_spatial = site_idx%VOLSPACE
 */
site_idx inline get_site_idx(const st_idx in)
{
  return in.space + VOLSPACE * in.time;
}
st_idx inline get_st_idx_from_site_idx(const site_idx in)
{
  st_idx tmp;
  tmp.space = in%VOLSPACE;
  tmp.time = in/VOLSPACE;
  return tmp;
}

/**
 * (st_idx, dir_idx) -> link_idx and
 * link_idx -> st_idx, respectively. 
 *using the convention:
 *link_idx = mu + NDIM * site_idx
 */
link_idx inline get_link_idx(const dir_idx mu, const st_idx in)
{
  return mu + NDIM * get_site_idx(in);
}
st_idx inline get_st_idx_from_link_idx(const link_idx in)
{
  st_idx tmp;
  site_idx idx_tmp = in/NDIM;
  tmp = get_st_idx_from_site_idx(idx_tmp);
  return tmp;
}

////////////////////////////////////////////////////////////////
// Accessing the lattice points
////////////////////////////////////////////////////////////////

/** returns all coordinates from a given site_idx */
coord_full get_coord_full(const site_idx in){
  st_idx tmp;
  tmp = get_st_idx_from_site_idx(in);
  coord_spatial tmp2;
  tmp2 = get_coord_spatial(tmp.space);
  coord_full res;
  res.w = tmp.time;
  res.x = tmp2.x;
  res.y = tmp2.y;
  res.z = tmp2.z;
  return res;
}

////////////////////////////////////////////////////////////////////
// Get positions of neighbours
// These functions do not rely explicitely on geometric convetions
////////////////////////////////////////////////////////////////////

/** returns neighbor in time direction given a temporal coordinate */
coord_temporal get_neighbor_temporal(const coord_temporal ntime){
  return (ntime + 1) % NTIME;
}

/** returns neighbor in time direction given a temporal coordinate */
coord_temporal get_lower_neighbor_temporal(const coord_temporal ntime){
  return (ntime - 1 + NTIME) % NTIME;
}

///////////////////////////////////////////////////////////////////
// Old functions
//////////////////////////////////////////////////////////////////

int inline get_global_pos(int spacepos, int t)
{
	return spacepos + VOLSPACE * t;
}

int inline get_global_link_pos(int mu, int spacepos, int t)
{
	return mu + NDIM * get_global_pos(spacepos, t);
}

//old version, this can be deleted soon:
int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t)
{
	return c + 2 * a + 2 * NC * b + 2 * NC * NC * mu + 2 * NC * NC * NDIM * spacepos + 2 * NC * NC * NDIM * VOLSPACE * t;
}

int inline ocl_su3matrix_element(int a, int b)
{
	return a + NC * b;
}

/**
 * An index type storing space (index) and time (coordinate) seperate.
 *
 * @todo Check name and types for suitability
 */
typedef struct {
	size_t space;
	uint time;
} st_index;

//it is assumed that idx iterates only over half the number of sites
st_index get_even_site(const size_t idx)
{
	uint4 tmp;
	tmp.x = idx;
	tmp.w = (int)(idx / (VOLSPACE / 2));
	tmp.x -= tmp.w * VOLSPACE / 2;
	tmp.z = (int)(tmp.x / (NSPACE * NSPACE / 2));
	tmp.x -= tmp.z * NSPACE * NSPACE / 2;
	tmp.y = (int)(tmp.x / NSPACE);
	tmp.x -= tmp.y * NSPACE;

	st_index res;
	res.space = (uint)((tmp.z + tmp.w) % 2) * (1 + 2 * tmp.x - (uint) (2 * tmp.x / NSPACE)) + (uint)((tmp.w + tmp.z + 1) % 2) * (2 * tmp.x + (uint) (2 * tmp.x / NSPACE)) + 2 * NSPACE * tmp.y + NSPACE * NSPACE * tmp.z;
	res.time = tmp.w;

	return res;
}

//it is assumed that idx iterates only over half the number of sites
st_index get_odd_site(const size_t idx)
{
	uint4 tmp;
	tmp.x = idx;
	tmp.w = (int)(idx / (VOLSPACE / 2));
	tmp.x -= tmp.w * VOLSPACE / 2;
	tmp.z = (int)(tmp.x / (NSPACE * NSPACE / 2));
	tmp.x -= tmp.z * NSPACE * NSPACE / 2;
	tmp.y = (int)(tmp.x / NSPACE);
	tmp.x -= tmp.y * NSPACE;

	st_index res;
	res.space =  (uint)((tmp.z + tmp.w + 1) % 2) * (1 + 2 * tmp.x - (uint) (2 * tmp.x / NSPACE)) + (uint)((tmp.w + tmp.z) % 2) * (2 * tmp.x + (uint) (2 * tmp.x / NSPACE)) + 2 * NSPACE * tmp.y + NSPACE * NSPACE * tmp.z;
	res.time = tmp.w;

	return res;
}

size_t get_nspace(const uint3 coord)
{
	return (coord.x +  NSPACE * coord.y + NSPACE * NSPACE * coord.z);
}

int get_spacecoord(const int nspace, const int dir)
{
	int res = nspace / (NSPACE * NSPACE);
	if(dir == 3) return res;
	int acc = res;
	res = nspace / NSPACE - NSPACE * acc;
	if(dir == 2) return res;
	acc = NSPACE * acc + res;
	res = nspace - NSPACE * acc;
	return res;
}

uint3 get_allspacecoord(const size_t nspace)
{
	uint3 coord;
	coord.z = nspace / NSPACE / NSPACE;
	int acc = coord.z;
	coord.y = nspace / NSPACE - NSPACE * acc;
	acc = NSPACE * acc + coord.y;
	coord.x = nspace - NSPACE * acc;
	return coord;
}

size_t get_neighbor(const size_t nspace, const uint dir)
{
	uint3 coord = get_allspacecoord(nspace);
	switch(dir) {
		case 1:
			coord.x = (coord.x + 1) % NSPACE;
			break;
		case 2:
			coord.y = (coord.y + 1) % NSPACE;
			break;
		case 3:
			coord.z = (coord.z + 1) % NSPACE;
			break;
	}
	return get_nspace(coord);
}

size_t get_lower_neighbor(const size_t nspace, const uint dir)
{
	uint3 coord = get_allspacecoord(nspace);
	switch(dir) {
		case 1:
			coord.x = (coord.x + NSPACE - 1) % NSPACE;
			break;
		case 2:
			coord.y = (coord.y + NSPACE - 1) % NSPACE;
			break;
		case 3:
			coord.z = (coord.z + NSPACE - 1) % NSPACE;
			break;
	}
	return get_nspace(coord);
}

int spinor_color(int spinor_element)
{
	return (int)(spinor_element / NSPIN);
}

int spinor_spin(int spinor_element, int color)
{
	return spinor_element - NSPIN * color;
}

int spinor_element(int alpha, int color)
{
	return alpha + NSPIN * color;
}

int get_n_eoprec(int spacepos, int timepos)
{
	return get_global_pos(spacepos, timepos) / 2;
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec)
{
	return alpha + NSPIN * color + NSPIN * NC * n_eoprec;
}

int spinor_field_element(int alpha, int color, int nspace, int t)
{
	return alpha + NSPIN * color + NSPIN * NC * (get_global_pos(nspace, t));
}
