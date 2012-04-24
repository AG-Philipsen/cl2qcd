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
 * #define XDIR 1
 * #define YDIR 2
 * #define ZDIR 3    (same)
 * A spatial idx is adressed as spatial_idx(x,y,z) = x * NS^(3-XDIR) + y * NS^(3-YDIR) + z * NS^(3-ZDIR) (different)
 * A site idx is addressed as site_idx(x,y,z,t) = spatial_idx(x,y,z) + t * VOLSPACE (same)
 * A link idx is addressed as an array with 2 indices.
 * @NOTE: This is more consistent then our convention because it yields
 *    site_idx = x * NS^(3-XDIR) + y * NS^(3-YDIR) + z * NS^(3-ZDIR) + t * NS^(3-TDIR)
 * whereas ours gives
 *    site_idx = x * NS^(XDIR-1) + y * NS^(YDIR-1) + z * NS^(ZDIR-1) + t * NS^(3-TDIR)
 */

////////////////////////////////////////////////////////////////
// Functions relying explicitely on the geometric conventions
// defined above (w/o even-odd functions!!)
////////////////////////////////////////////////////////////////

/**
 * with this set to false or true, one can switch between our original convention and
 * the one from tmlqcd.
 * our original:
 * spatial_idx = x + NS * y + NS*NS * z
 * tmlqcd:
 * spatial_idx = z + NS * y + NS*NS * x
 * NOTE: the ifs and elses used here should be removed by the compiler
 *       Nevertheless, one could also change to a permanent convention here
 */
#define TMLQCD_CONV false

/** spatial coordinates <-> spatial_idx
 *@todo this can be generalize using the definitions of the spatial directions
 *see  http://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
 */
spatial_idx get_spatial_idx(const coord_spatial coord)
{
	bool tmp = TMLQCD_CONV;
	if( tmp ) {
		return (coord.z +  NSPACE * coord.y + NSPACE * NSPACE * coord.x);
	} else {
		return (coord.x +  NSPACE * coord.y + NSPACE * NSPACE * coord.z);
	}
}
coord_spatial get_coord_spatial(const spatial_idx nspace)
{
	coord_spatial coord;
	bool tmp = TMLQCD_CONV;
	if( tmp ) {
		coord.x = nspace / NSPACE / NSPACE;
		uint acc = coord.x;
		coord.y = nspace / NSPACE - NSPACE * acc;
		acc = NSPACE * acc + coord.y;
		coord.z = nspace - NSPACE * acc;
	} else {
		coord.z = nspace / NSPACE / NSPACE;
		uint acc = coord.z;
		coord.y = nspace / NSPACE - NSPACE * acc;
		acc = NSPACE * acc + coord.y;
		coord.x = nspace - NSPACE * acc;
	}
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
	tmp.space = in % VOLSPACE;
	tmp.time = in / VOLSPACE;
	return tmp;
}

/** returns eo-vector component from st_idx */
site_idx get_eo_site_idx_from_st_idx(st_idx in);
st_idx get_odd_st_idx(const site_idx idx);
st_idx get_even_st_idx(const site_idx idx);

/**
 * (st_idx, dir_idx) -> link_idx and
 * link_idx -> st_idx, respectively.
 *using the convention:
 *link_idx = mu + NDIM * site_idx
 */
link_idx inline get_link_idx_SOA(const dir_idx mu, const st_idx in)
{
	const uint3 space = get_coord_spatial(in.space);
	// check if the site is odd (either spacepos or t odd)
	// if yes offset everything by half the number of sites (number of sites is always even...)
	site_idx odd = (space.x ^ space.y ^ space.z ^ in.time) & 0x1 ? (VOL4D / 2) : 0;
	return mu * VOL4D + odd + get_eo_site_idx_from_st_idx(in);
}
link_idx inline get_link_idx_AOS(const dir_idx mu, const st_idx in)
{
	return mu + NDIM * get_site_idx(in);
}
link_idx inline get_link_idx(const dir_idx mu, const st_idx in)
{
#ifdef _USE_SOA_
	return get_link_idx_SOA(mu, in);
#else
	return get_link_idx_AOS(mu, in);
#endif
}

dir_idx inline get_dir_idx_from_link_idx(const link_idx in)
{
#ifdef _USE_SOA_
	return in / NDIM;
#else
	return in % NDIM;
#endif
}

////////////////////////////////////////////////////////////////
// Accessing the lattice points
////////////////////////////////////////////////////////////////

/** returns all coordinates from a given site_idx */
coord_full get_coord_full(const site_idx in)
{
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

////////////////////////////////////////////////////////////////
// Even/Odd functions
// NOTE: These explicitely rely on the geometric conventions
// and should be reconsidered if the former are changed!!
////////////////////////////////////////////////////////////////

/**
 * Even-Odd preconditioning means applying a chessboard to the lattice,
 * and distinguish between black and white (even/odd) sites.
 * this can be done in a 2D plane. Going to 4D then means to repeat this
 * 2D plane alternating even and odd sites.
 *
 * Schematically one has (assume 4x4 plane, x: even, o: odd):
 * oxox
 * xoxo
 * oxox
 * xoxo
 *
 * Take a look at the coordinates (call them w (ordinate) and t (abscise) ):
 * w/t | 0 | 1 | 2 | 3 |
 *   0 | x | o | x | o |
 *   1 | o | x | o | x |
 *   2 | x | o | x | o |
 *   3 | o | x | o | x |
 *
 * Assume that a coordinate in the plane is calculated as pos = t + w*4
 * Then, the coordinates of the even and odd sites are:
 * even|odd
 * 0|1
 * 2|3
 * 5|4
 * 7|6
 * 9|8
 * 11|10
 * 12|13
 * 14|15
 *
 * One can see a regular pattern, which can be described by the functions
 * f_even:2*w + int(2w/4)
 * f_odd: 1 + 2*w - int(2w/4)
 * Thus, given an index pair (w,t) one can map it to an even or odd site
 * by pos = f_even/odd(w) + 4*2*t, assuming that t runs only from 0 to 2
 * (half lattice size).
 * Extending this to 4D (call additional coordinates g1,g2) can be done
 * by alternating between f_even and f_odd, depending on g1 and g2.
 * E.g. one can do that via (g1 + g2)%2 and (g1 + g2 + 1)%2 prefactors.
 * Extending these consideratios to arbitrary lattice extend is straightforward.
 * NOTE: In general, a 4D site is even or odd if (x+y+z+t)%2 is 0 or 1.
 * However, in general one wants to have a loop over all even or odd sites.
 * Then this condition is not applicable, if not one loops over ALL sites and
 * just leaves out the ones where the condition is not fulfilled.
 * Since we do not want to create a possibly big table to map loop-variable
 * -> even/odd site we will use the pattern above for that!
 */

/** returns eo-vector component from st_idx */
site_idx get_eo_site_idx_from_st_idx(st_idx in)
{
	return get_site_idx(in) / 2;
}

/** the functions mentioned above
 * @todo this may be generalized because it relys on the current geometric conventions!!
 */
site_idx calc_even_spatial_idx(coord_full in)
{
	bool switcher = TMLQCD_CONV;
	if(switcher) {
		return
		  (uint)((in.x + in.w    ) % 2) * (1 + 2 * in.z - (uint) (2 * in.z / NSPACE)) +
		  (uint)((in.w + in.x + 1) % 2) * (    2 * in.z + (uint) (2 * in.z / NSPACE)) +
		  2 * NSPACE * in.y +
		  NSPACE * NSPACE * in.x;
	} else {
		return
		  (uint)((in.z + in.w    ) % 2) * (1 + 2 * in.x - (uint) (2 * in.x / NSPACE)) +
		  (uint)((in.w + in.z + 1) % 2) * (    2 * in.x + (uint) (2 * in.x / NSPACE)) +
		  2 * NSPACE * in.y +
		  NSPACE * NSPACE * in.z;
	}
}
site_idx calc_odd_spatial_idx(coord_full in)
{
	bool switcher = TMLQCD_CONV;
	if(switcher) {
		return
		  (uint)((in.x + in.w + 1) % 2) * (1 + 2 * in.z - (uint) (2 * in.z / NSPACE)) +
		  (uint)((in.w + in.x    ) % 2) * (    2 * in.z + (uint) (2 * in.z / NSPACE)) +
		  2 * NSPACE * in.y +
		  NSPACE * NSPACE * in.x;
	} else {
		return
		  (uint)((in.z + in.w + 1) % 2) * (1 + 2 * in.x - (uint) (2 * in.x / NSPACE)) +
		  (uint)((in.w + in.z    ) % 2) * (    2 * in.x + (uint) (2 * in.x / NSPACE)) +
		  2 * NSPACE * in.y +
		  NSPACE * NSPACE * in.z;
	}
}

/**
 * this takes a eo_site_idx (0..VOL4D/2) and returns its 4 coordinates
 * under the assumption that even-odd preconditioning is applied in the
 * x-y-plane as described above.
 * This is moved to the z-y plane if the tmlqcd conventions are used.
 * This is done for convenience since x and y direction have
 * the same extent.
 * Use the convention:
 *site_idx = x + y*NS + z*NS*NS + t*NS*NS*NS
 *= site_idx_spatial + t*VOLSPACE
 *and
 *spatial_idx = x * NS^XDIR-1 + y * NS^YDIR-1 + z * NS^ZDIR-1
 *= x + y * NS + z * NS^2
 *
 * Then one can "dissect" the site_idx i according to
 * t= i/VOLSPACE/2
 * z = (i-t*VOLSPACE/2)/(NS*NS/2)
 * y = (i-t*VOLSPACE/2 - z*NS*NS/2) / NS
 * x = (i-t*VOLSPACE/2 - z*NS*NS/2 - y*NS)
 * As mentioned above, y is taken to run from 0..NS/2
 */
coord_full dissect_eo_site_idx(const site_idx idx)
{
	coord_full tmp;
	bool switcher = TMLQCD_CONV;
	if(switcher) {
		tmp.z = idx;
		tmp.w = (int)(idx / (VOLSPACE / 2));
		tmp.z -= tmp.w * VOLSPACE / 2;
		tmp.x = (int)(tmp.z / (NSPACE * NSPACE / 2));
		tmp.z -= tmp.x * NSPACE * NSPACE / 2;
		tmp.y = (int)(tmp.z / NSPACE);
		tmp.z -= tmp.y * NSPACE;
	} else {
		tmp.x = idx;
		tmp.w = (int)(idx / (VOLSPACE / 2));
		tmp.x -= tmp.w * VOLSPACE / 2;
		tmp.z = (int)(tmp.x / (NSPACE * NSPACE / 2));
		tmp.x -= tmp.z * NSPACE * NSPACE / 2;
		tmp.y = (int)(tmp.x / NSPACE);
		tmp.x -= tmp.y * NSPACE;
	}
	return tmp;
}

/** given an eo site_idx (0..VOL4D/2), returns corresponding even site_idx */
st_idx get_even_st_idx(const site_idx idx)
{
	coord_full tmp = dissect_eo_site_idx(idx);
	st_idx res;
	res.space = calc_even_spatial_idx(tmp);
	res.time = tmp.w;
	return res;
}

/** given an eo site_idx (0..VOL4D/2), returns corresponding odd site_idx */
st_idx get_odd_st_idx(const site_idx idx)
{
	coord_full tmp = dissect_eo_site_idx(idx);
	st_idx res;
	res.space = calc_odd_spatial_idx(tmp);
	res.time = tmp.w;
	return res;
}

////////////////////////////////////////////////////////////////////
// Get positions of neighbours
// These functions do not rely explicitely on geometric convetions
////////////////////////////////////////////////////////////////////

/** returns neighbor in time direction given a temporal coordinate */
coord_temporal get_neighbor_temporal(const coord_temporal ntime)
{
	return (ntime + 1) % NTIME;
}

/** returns neighbor in time direction given a temporal coordinate */
coord_temporal get_lower_neighbor_temporal(const coord_temporal ntime)
{
	return (ntime - 1 + NTIME) % NTIME;
}

/** returns idx of neighbor in spatial direction dir given a spatial idx */
site_idx get_neighbor_spatial(const spatial_idx nspace, const dir_idx dir)
{
	coord_spatial coord = get_coord_spatial(nspace);
	switch(dir) {
		case XDIR:
			coord.x = (coord.x + 1) % NSPACE;
			break;
		case YDIR:
			coord.y = (coord.y + 1) % NSPACE;
			break;
		case ZDIR:
			coord.z = (coord.z + 1) % NSPACE;
			break;
	}
	return get_spatial_idx(coord);
}

/** returns idx of lower neighbor in spatial direction dir given a spatial idx */
site_idx get_lower_neighbor_spatial(const spatial_idx nspace, const dir_idx dir)
{
	coord_spatial coord = get_coord_spatial(nspace);
	switch(dir) {
		case XDIR:
			coord.x = (coord.x - 1 + NSPACE) % NSPACE;
			break;
		case YDIR:
			coord.y = (coord.y - 1 + NSPACE) % NSPACE;
			break;
		case ZDIR:
			coord.z = (coord.z - 1 + NSPACE) % NSPACE;
			break;
	}
	return get_spatial_idx(coord);
}

/** returns the st_idx of the neighbor in direction dir given a st_idx. */
st_idx get_neighbor_from_st_idx(const st_idx in, const dir_idx dir)
{
	st_idx tmp = in;
	if(dir == TDIR) tmp.time = get_neighbor_temporal(in.time);
	else tmp.space = get_neighbor_spatial(in.space, dir);
	return tmp;
}

/** returns the st_idx of the lower neighbor in direction dir given a st_idx. */
st_idx get_lower_neighbor_from_st_idx(const st_idx in, const dir_idx dir)
{
	st_idx tmp = in;
	if(dir == TDIR) tmp.time = get_lower_neighbor_temporal(in.time);
	else tmp.space = get_lower_neighbor_spatial(in.space, dir);
	return tmp;
}

////////////////////////////////////////////////////////////////
// Accessing spinor elements
////////////////////////////////////////////////////////////////

uint spinor_color(uint spinor_element)
{
	return (uint)(spinor_element / NSPIN);
}

uint spinor_spin(uint spinor_element, uint color)
{
	return spinor_element - NSPIN * color;
}

uint spinor_element(uint alpha, uint color)
{
	return alpha + NSPIN * color;
}

///////////////////////////////////////////////////////////////////
// Older functions and types which still appear in the code
// CP: I simply directed them back to the newer functions
//////////////////////////////////////////////////////////////////

int inline get_global_pos(int spacepos, int t)
{
	st_idx tmp = {spacepos, t};
	return get_site_idx(tmp);
}

int inline get_global_link_pos(int mu, int spacepos, int t)
{
	st_idx tmp = {spacepos, t};
	return get_link_idx(mu, tmp);
}

/**
 * An index type storing space (index) and time (coordinate) seperate.
 *
 * @NOTE for convenience I changed this simply to "st_idx"
 */
typedef st_idx st_index;

st_index inline get_even_site(const size_t idx)
{
	return get_even_st_idx(idx);
}

st_index inline get_odd_site(const size_t idx)
{
	return get_odd_st_idx(idx);
}

size_t inline get_nspace(const uint3 coord)
{
	return (size_t) get_spatial_idx(coord);
}

uint3  inline get_allspacecoord(const size_t nspace)
{
	return get_coord_spatial(nspace);
}

size_t inline get_neighbor(const size_t nspace, const uint dir)
{
	return get_neighbor_spatial(nspace, dir);
}

size_t inline get_lower_neighbor(const size_t nspace, const uint dir)
{
	return get_lower_neighbor_spatial(nspace, dir);
}

int inline get_n_eoprec(int spacepos, int timepos)
{
	st_idx tmp = {spacepos, timepos};
	return get_eo_site_idx_from_st_idx(tmp);
}

