__kernel void prng_ranlux_init(uint ins, __global rngStateStorageType * const restrict states)
{
	ranluxcl_initialization(ins, states);
}
