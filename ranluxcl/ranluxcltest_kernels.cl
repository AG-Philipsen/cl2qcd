/**
 * Copyright (c) 2013 Matthias Bach
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
 *
 * NOTE: The code contained in this file was developed by external developers
 *       and the copyright and license statements above refer to the work
 *       that was done to include the third party code into CL2QCD.
 */

// If the include does not work then the contents of the ranluxcl.cl file can just be
// pasted here at the top instead.
#include "ranluxcl.cl"

__kernel void Kernel_Ranluxcl_Init(private int ins, global float4* ranluxcltab)
{
    ranluxcl_initialization(ins, ranluxcltab);
}

__kernel void Kernel_PRN(private int KernelCycles, global float4* ranluxcltab, global float* PRNs)
{
    // Downloading ranluxcltab. The state of RANLUXCL is stored in ranluxclstate.
    ranluxcl_state_t ranluxclstate;
    ranluxcl_download_seed(&ranluxclstate, ranluxcltab);

    float4 randomnr;

    // Generate some numbers
    for (int i = 0; i < KernelCycles; i += 4)
        randomnr = ranluxcl(&ranluxclstate);

    // Uploading only last number generated.
    PRNs[get_global_id(0)] = randomnr.w;

    // Uploading ranluxcltab
    ranluxcl_upload_seed(&ranluxclstate, ranluxcltab);
}
