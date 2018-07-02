/*
 * Copyright (c) 2011-2013 Matthias Bach
 * Copyright (c) 2011 Christopher Pinke
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file stout-smearing-fermion-force
 */

// this is done after the tmlqcd-analogue (stout_smear_force.c)
// I did not check this yet, but I think since the calculation of the stout-smeared force is iterative, one should most
// likely split this into 2 kernels...

// this kernel is meant to be the one that updates the force with all needed ingredients already computed
//    it takes the force field (out) and <fill in>
__kernel void stout_smear_fermion_force(__global aeStorageType* const restrict out)
{
    int local_size  = get_local_size(0);
    int global_size = get_global_size(0);
    int id          = get_global_id(0);
    int loc_idx     = get_local_id(0);
    int num_groups  = get_num_groups(0);
    int group_id    = get_group_id(0);
}
