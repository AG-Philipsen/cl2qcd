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

//dispensable
//site = pos + VOLSPACE*t = x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
/*
int inline ocl_gaugefield_element(int mu, int spacepos, int t)
{
return get_global_link_pos(mu, spacepos, t);
}
*/

//old version, this can be deleted soon:
int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t)
{
#ifdef _RECONSTRUCT_TWELVE_
	return c + 2 * a + 2 * (NC - 1) * b + 2 * NC * (NC - 1) * mu + 2 * NC * (NC - 1) * NDIM * spacepos + 2 * NC * (NC - 1) * NDIM * VOLSPACE * t;
#else
	return c + 2 * a + 2 * NC * b + 2 * NC * NC * mu + 2 * NC * NC * NDIM * spacepos + 2 * NC * NC * NDIM * VOLSPACE * t;
#endif
}

int inline ocl_su3matrix_element(int a, int b)
{
#ifdef _RECONSTRUCT_TWELVE_
	return a + (NC - 1) * b;
#else
	return a + NC * b;
#endif
}


//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
// int inline ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t)
// {
// #ifdef _RECONSTRUCT_TWELVE_
//  return c + 2*a + 2*(NC-1)*b+2*NC*(NC-1)*mu+2*NC*(NC-1)*NDIM*spacepos+2*NC*(NC-1)*NDIM*VOLSPACE*t;
// #else
//  return c + 2*a + 2*NC*b+2*NC*NC*mu+2*NC*NC*NDIM*spacepos+2*NC*NC*NDIM*VOLSPACE*t;
// #endif
// }

//dispensable
/*
int inline ocl_su3matrix_element(int a, int b)
{
#ifdef _RECONSTRUCT_TWELVE_
  return a + (NC-1)*b;
#else
  return a + NC*b;
#endif
}
*/

//it is assumed that idx iterates only over half the number of sites
void inline get_even_site(int idx, int * out_space, int * out_t)
{
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;
	(*out_space) =  (int)((z + t) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z + 1) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites
void inline get_odd_site(int idx, int * out_space, int * out_t)
{
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;
	(*out_space) =  (int)((z + t + 1) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}

int get_nspace(const int* coord)
{
	int n = coord[1] +  NSPACE * coord[2] + NSPACE * NSPACE * coord[3];
	return n;
}

int get_spacecoord(const int nspace, const int dir)
{
	int res = convert_int(nspace / (NSPACE * NSPACE));
	if(dir == 3) return res;
	int acc = res;
	res = convert_int(nspace / (NSPACE)) - NSPACE * acc;
	if(dir == 2) return res;
	acc = NSPACE * acc + res;
	res = nspace - NSPACE * acc;
	return res;
}

void get_allspacecoord(const int nspace, int coord[NDIM])
{
	int res = convert_int(nspace / (NSPACE * NSPACE));
	coord[3] = res;
	int acc = res;
	res = convert_int(nspace / (NSPACE)) - NSPACE * acc;
	coord[2] = res;
	acc = NSPACE * acc + res;
	res = nspace - NSPACE * acc;
	coord[1] = res;
}

int get_neighbor(const int nspace, const int dir)
{
	int coord[NDIM];
	get_allspacecoord(nspace, coord);
	coord[dir] = (coord[dir] + 1) % NSPACE;
	return get_nspace(coord);
}

int get_lower_neighbor(const int nspace, int const dir)
{
	int coord[NDIM];
	get_allspacecoord(nspace, coord);
	coord[dir] = (coord[dir] - 1 + NSPACE) % NSPACE;
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
	return (int)((get_global_pos(spacepos, timepos)) / 2);
}

int eoprec_spinor_field_element(int alpha, int color, int n_eoprec)
{
	return alpha + NSPIN * color + NSPIN * NC * n_eoprec;
}

int spinor_field_element(int alpha, int color, int nspace, int t)
{
	return alpha + NSPIN * color + NSPIN * NC * (get_global_pos(nspace, t));
}
