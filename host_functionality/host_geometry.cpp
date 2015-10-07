/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "host_geometry.h"

#include "../meta/util.hpp"

int get_nspace(int* coord, const meta::Inputparameters& params)
{
	return get_nspace(coord, params.get_ntime(), params.get_nspace() );
}

int get_nspace(int* coord, const int nt, const int ns)
{
	int n = 0;
	n = ns * ns * coord[3] + ns * coord[2] + coord[1];
	//for(int j = 1; j < NDIM; j++) n += pow(params->get_ns(), j - 1) * coord[j];
	return n;
}

//functions that have explicit spatial and temporal positions in them
//make it:
//site = pos + VOLSPACE*t =  x + y*NSPACE + z*NSPACE*NSPACE + VOLSPACE*t
//idx = mu + NDIM*site
int get_global_pos(int spacepos, int t, const meta::Inputparameters& params)
{
	return spacepos + meta::get_volspace(params) * t;
}
int get_global_pos(int spacepos, int t, const int nt, const int ns)
{
	return spacepos + meta::get_volspace(ns) * t;
}


int get_global_pos(size_4 cart, const meta::Inputparameters& params)
{
    int space_coords[4];
    space_coords[0] = 0;
    space_coords[1] = cart.x;
    space_coords[2] = cart.y;
    space_coords[3] = cart.z;
    auto const spatial = get_nspace(space_coords, params);
    return get_global_pos(spatial, cart.t, params);
}

int get_global_pos(size_4 cart, const int nt, const int ns)
{
    int space_coords[4];
    space_coords[0] = 0;
    space_coords[1] = cart.x;
    space_coords[2] = cart.y;
    space_coords[3] = cart.z;
    auto const spatial = get_nspace(space_coords, nt, ns);
    return get_global_pos(spatial, cart.t, nt, ns);
}

int get_global_link_pos(int mu, int spacepos, int t, const meta::Inputparameters& params)
{
	return mu + NDIM * get_global_pos(spacepos, t, params);
}
int get_global_link_pos(int mu, int spacepos, int t, const int nt, const int ns)
{
	return mu + NDIM * get_global_pos(spacepos, t, nt, ns);
}
int get_global_link_pos(int mu, size_4 cart, const meta::Inputparameters& params)
{
	int space_coords[4];
	space_coords[0] = 0;
	space_coords[1] = cart.x;
	space_coords[2] = cart.y;
	space_coords[3] = cart.z;
	auto const spatial = get_nspace(space_coords, params);
	return get_global_link_pos(mu, spatial, cart.t, params);
}
int get_global_link_pos(int mu, size_4 cart, const int nt, const int ns)
{
	int space_coords[4];
	space_coords[0] = 0;
	space_coords[1] = cart.x;
	space_coords[2] = cart.y;
	space_coords[3] = cart.z;
	auto const spatial = get_nspace(space_coords, nt, ns);
	return get_global_link_pos(mu, spatial, cart.t, nt, ns);
}

size_t get_su3_idx_ildg_format(size_t n, size_t m, size_t x, size_t y, size_t z, size_t t, size_t mu, const int nt, const int ns)
{
	//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
	size_t link_idx_ildg = get_link_idx_ildg_format(x, y, z, t, mu, nt, ns);
	//the 2 is introduced because we are dealing with complex numbers here..
	return 2 * n + 2 * m * NC + link_idx_ildg * NC * NC * 2;
}

size_t get_link_idx_ildg_format(size_t x, size_t y, size_t z, size_t t, size_t mu, const int nt, const int ns)
{
	const size_t VOLSPACE = meta::get_volspace(ns);
	size_t spacepos_ildg = z + y * ns + x * ns * ns;
	size_t link_idx_ildg = mu + spacepos_ildg * NDIM + t * VOLSPACE * NDIM;
	return link_idx_ildg;
}

int get_source_pos_spatial(const meta::Inputparameters& params)
{
	int coord [4];
	coord[1] = params.get_source_x();
	coord[2] = params.get_source_y();
	coord[3] = params.get_source_z();
	return get_nspace(coord, params);
}

