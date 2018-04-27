/*
 * Copyright (c) 2011 Christopher Pinke
 * Copyright (c) 2012,2013 Matthias Bach
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

/** @file
 * Device code for the heatbath overrelaxation
 */

__kernel void overrelax_odd(__global Matrixsu3StorageType* const restrict gaugefield, const int mu,
                            __global rngStateStorageType* const restrict rngStates)
{
    prng_state rnd;
    prng_loadState(&rnd, rngStates);

    PARALLEL_FOR (id, VOL4D_LOCAL / 2) {
        st_index pos = get_odd_st_idx_local(id);
        perform_overrelaxing(gaugefield, mu, &rnd, pos.space, pos.time);
    }

    prng_storeState(rngStates, &rnd);
}
