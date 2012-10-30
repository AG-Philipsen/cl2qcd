/**
 @file fermion-observables: chiral condensate
*/

//this is pbp in the standard version
__kernel void pbp_std(__global const spinor * const restrict phi,__global const spinor * const restrict b,  __global hmc_float * const restrict out)
{
    int local_size = get_local_size(0);
    int global_size = get_global_size(0);
    int id = get_global_id(0);
    int loc_idx = get_local_id(0);
    int num_groups = get_num_groups(0);
    int group_id = get_group_id (0);

    return;
}

//this is pbp in the one-end-trick version for twisted mass
__kernel void pbp_tm_one_end(__global const spinor * const restrict phi, __global const spinor * const restrict b,  __global hmc_float * const restrict out)
{
    int local_size = get_local_size(0);
    int global_size = get_global_size(0);
    int id = get_global_id(0);
    int loc_idx = get_local_id(0);
    int num_groups = get_num_groups(0);
    int group_id = get_group_id (0);

    return;
}
