/** @file
 * Device code for lattice geometry handling
 */

//opencl_geometry.cl

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
