/** @file
 * Kernel for bandwidth microbenchmark
 */

/**
 * Calculate the proper stride for SOA storage.
 *
 * \param The number of elements in the array.
 * \param The size of the datatype used for the array.
 * \return The proper stride in elements of the storage array.
 */
ulong calculateStride(const ulong elems, const ulong baseTypeSize)
{
	// Align stride to (N * 16 + 8) KiB
	// TODO this is optimal for AMD HD 5870, also adjust for others
	const ulong stride_bytes = ((elems * baseTypeSize + 0x1FFF) & 0xFFFFFFFFFFFFC000L) | 0x2000;
	const ulong stride_elems = stride_bytes / baseTypeSize;
	return stride_elems;
}

__kernel void copyFloat( __global hmc_float * const restrict out, __global const hmc_float * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		hmc_float tmp = in[i];
		out[i] = tmp;
	}
}

__kernel void copySU3( __global Matrixsu3 * const restrict out, __global const Matrixsu3 * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		Matrixsu3 tmp = in[i];
		out[i] = tmp;
	}
}

Matrixsu3 peekSU3SOA(__global const hmc_float * const restrict in, const size_t idx, const ulong stride)
{
	return (Matrixsu3) {
		{
			in[ 0 * stride + idx], in[ 1 * stride + idx]
		},
		{ in[ 2 * stride + idx], in[ 3 * stride + idx] },
		{ in[ 4 * stride + idx], in[ 5 * stride + idx] },
		{ in[ 6 * stride + idx], in[ 7 * stride + idx] },
		{ in[ 8 * stride + idx], in[ 9 * stride + idx] },
		{ in[10 * stride + idx], in[11 * stride + idx] },
		{ in[12 * stride + idx], in[13 * stride + idx] },
		{ in[14 * stride + idx], in[15 * stride + idx] },
		{ in[16 * stride + idx], in[17 * stride + idx] }
	};
}

void pokeSU3SOA(__global hmc_float * const restrict out, const size_t idx, const Matrixsu3 val, const ulong stride)
{
	out[ 0 * stride + idx] = val.e00.re;
	out[ 1 * stride + idx] = val.e00.im;
	out[ 2 * stride + idx] = val.e01.re;
	out[ 3 * stride + idx] = val.e01.im;
	out[ 4 * stride + idx] = val.e02.re;
	out[ 5 * stride + idx] = val.e02.im;
	out[ 6 * stride + idx] = val.e10.re;
	out[ 7 * stride + idx] = val.e10.im;
	out[ 8 * stride + idx] = val.e11.re;
	out[ 9 * stride + idx] = val.e11.im;
	out[10 * stride + idx] = val.e12.re;
	out[11 * stride + idx] = val.e12.im;
	out[12 * stride + idx] = val.e20.re;
	out[13 * stride + idx] = val.e20.im;
	out[14 * stride + idx] = val.e21.re;
	out[15 * stride + idx] = val.e21.im;
	out[16 * stride + idx] = val.e22.re;
	out[17 * stride + idx] = val.e22.im;
}

__kernel void copySU3SOA( __global hmc_float * const restrict out, __global const hmc_float * const restrict in, const ulong elems, const ulong threads_per_group )
{
	const size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	const size_t id = get_group_id(0) * threads_per_group + local_id;
	const ulong stride = calculateStride(elems, sizeof(hmc_float));

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		Matrixsu3 tmp = peekSU3SOA(in, i, stride);
		pokeSU3SOA(out, i, tmp, stride);
	}
}

Matrixsu3 peekSU3SOcplxA(__global const hmc_complex * const restrict in, const size_t idx, const ulong stride)
{
	return (Matrixsu3) {
		in[0 * stride + idx],
		in[1 * stride + idx],
		in[2 * stride + idx],
		in[3 * stride + idx],
		in[4 * stride + idx],
		in[5 * stride + idx],
		in[7 * stride + idx],
		in[8 * stride + idx],
		in[9 * stride + idx]
	};
}

void pokeSU3SOcplxA(__global hmc_complex * const restrict out, const size_t idx, const Matrixsu3 val, const ulong stride)
{
	out[0 * stride + idx] = val.e00;
	out[1 * stride + idx] = val.e01;
	out[2 * stride + idx] = val.e02;
	out[3 * stride + idx] = val.e10;
	out[4 * stride + idx] = val.e11;
	out[5 * stride + idx] = val.e12;
	out[6 * stride + idx] = val.e20;
	out[7 * stride + idx] = val.e21;
	out[8 * stride + idx] = val.e22;
}

__kernel void copySU3SOcplxA( __global hmc_complex * const restrict out, __global const hmc_complex * const restrict in, const ulong elems, const ulong threads_per_group )
{
	const size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	const size_t id = get_group_id(0) * threads_per_group + local_id;
	const ulong stride = calculateStride(elems, sizeof(hmc_complex));

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		Matrixsu3 tmp = peekSU3SOcplxA(in, i, stride);
		pokeSU3SOcplxA(out, i, tmp, stride);
	}
}

__kernel void copySpinor( __global spinor * const restrict out, __global const spinor * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		spinor tmp = in[i];
		out[i] = tmp;
	}
}

hmc_complex make_complex(const hmc_float left, const hmc_float right)
{
	return (hmc_complex) {
		left, right
	};
}

su3vec make_su3vec(const hmc_complex e0, const hmc_complex e1, const hmc_complex e2)
{
	return (su3vec) {
		e0, e1, e2
	};
}

spinor make_spinor(const su3vec e0, const su3vec e1, const su3vec e2, const su3vec e3)
{
	return (spinor) {
		e0, e1, e2, e3
	};
}

spinor peekSpinorSOA(__global const hmc_float * const restrict in, const size_t idx, const ulong stride)
{
	return make_spinor(
	  make_su3vec(
	    make_complex(in[ 0 * stride + idx], in[ 1 * stride + idx]),
	    make_complex(in[ 2 * stride + idx], in[ 3 * stride + idx]),
	    make_complex(in[ 4 * stride + idx], in[ 5 * stride + idx])
	  ), make_su3vec(
	    make_complex(in[ 6 * stride + idx], in[ 7 * stride + idx]),
	    make_complex(in[ 8 * stride + idx], in[ 9 * stride + idx]),
	    make_complex(in[10 * stride + idx], in[11 * stride + idx])
	  ), make_su3vec(
	    make_complex(in[12 * stride + idx], in[13 * stride + idx]),
	    make_complex(in[14 * stride + idx], in[15 * stride + idx]),
	    make_complex(in[16 * stride + idx], in[17 * stride + idx])
	  ), make_su3vec(
	    make_complex(in[18 * stride + idx], in[19 * stride + idx]),
	    make_complex(in[20 * stride + idx], in[21 * stride + idx]),
	    make_complex(in[22 * stride + idx], in[23 * stride + idx])
	  )
	);
}

void pokeSpinorSOA(__global hmc_float * const restrict out, const size_t idx, const spinor val, const ulong stride)
{
	// su3vec = 3 * cplx
	out[ 0 * stride + idx] = val.e0.e0.re;
	out[ 1 * stride + idx] = val.e0.e0.im;
	out[ 2 * stride + idx] = val.e0.e1.re;
	out[ 3 * stride + idx] = val.e0.e1.im;
	out[ 4 * stride + idx] = val.e0.e2.re;
	out[ 5 * stride + idx] = val.e0.e2.im;

	// su3vec = 3 * cplx
	out[ 6 * stride + idx] = val.e1.e0.re;
	out[ 7 * stride + idx] = val.e1.e0.im;
	out[ 8 * stride + idx] = val.e1.e1.re;
	out[ 9 * stride + idx] = val.e1.e1.im;
	out[10 * stride + idx] = val.e1.e2.re;
	out[11 * stride + idx] = val.e1.e2.im;

	// su3vec = 3 * cplx
	out[12 * stride + idx] = val.e2.e0.re;
	out[13 * stride + idx] = val.e2.e0.im;
	out[14 * stride + idx] = val.e2.e1.re;
	out[15 * stride + idx] = val.e2.e1.im;
	out[16 * stride + idx] = val.e2.e2.re;
	out[17 * stride + idx] = val.e2.e2.im;

	// su3vec = 3 * cplx
	out[18 * stride + idx] = val.e3.e0.re;
	out[19 * stride + idx] = val.e3.e0.im;
	out[20 * stride + idx] = val.e3.e1.re;
	out[21 * stride + idx] = val.e3.e1.im;
	out[22 * stride + idx] = val.e3.e2.re;
	out[23 * stride + idx] = val.e3.e2.im;
}

__kernel void copySpinorSOA( __global hmc_float * const restrict out, __global const hmc_float * const restrict in, const ulong elems, const ulong threads_per_group )
{
	const size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	size_t id = get_group_id(0) * threads_per_group + local_id;
	const ulong stride = calculateStride(elems, sizeof(hmc_float));
	for( size_t i = id; i < elems;  i += get_global_size(0) ) {
		spinor tmp = peekSpinorSOA(in, i, stride);
		pokeSpinorSOA(out, i, tmp, stride);
	}
}

spinor peekSpinorSOcplxA(__global const hmc_complex * const restrict in, const size_t idx, const ulong stride)
{
	return (spinor) {
		{
			// su3vec = 3 * cplx
			in[ 0 * stride + idx],
			in[ 1 * stride + idx],
			in[ 2 * stride + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 3 * stride + idx],
			in[ 4 * stride + idx],
			in[ 5 * stride + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 6 * stride + idx],
			in[ 7 * stride + idx],
			in[ 8 * stride + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 9 * stride + idx],
			in[10 * stride + idx],
			in[11 * stride + idx]
		}
	};
}

void pokeSpinorSOcplxA(__global hmc_complex * const restrict out, const size_t idx, const spinor val, const ulong stride)
{
	// su3vec = 3 * cplx
	out[ 0 * stride + idx] = val.e0.e0;
	out[ 1 * stride + idx] = val.e0.e1;
	out[ 2 * stride + idx] = val.e0.e2;

	// su3vec = 3 * cplx
	out[ 3 * stride + idx] = val.e1.e0;
	out[ 4 * stride + idx] = val.e1.e1;
	out[ 5 * stride + idx] = val.e1.e2;

	// su3vec = 3 * cplx
	out[ 6 * stride + idx] = val.e2.e0;
	out[ 7 * stride + idx] = val.e2.e1;
	out[ 8 * stride + idx] = val.e2.e2;

	// su3vec = 3 * cplx
	out[ 9 * stride + idx] = val.e3.e0;
	out[10 * stride + idx] = val.e3.e1;
	out[11 * stride + idx] = val.e3.e2;
}

__kernel void copySpinorSOcplxA( __global hmc_complex * const restrict out, __global const hmc_complex * const restrict in, const ulong elems, const ulong threads_per_group )
{
	size_t local_id = get_local_id(0);

	if( local_id >= threads_per_group )
		return;

	const size_t id = get_group_id(0) * threads_per_group + local_id;
	const ulong stride = calculateStride(elems, sizeof(hmc_complex));

	for( size_t i = id; i < elems;  i += threads_per_group * get_num_groups(0) ) {
		spinor tmp = peekSpinorSOcplxA(in, i, stride);
		pokeSpinorSOcplxA(out, i, tmp, stride);
	}
}

// FIXME for now simply assume we have 128 threads
__attribute__((reqd_work_group_size(128, 1, 1)))
__kernel void copySpinorLocal( __global spinor * const restrict out, __global const spinor * const restrict in, const ulong elems, const ulong threads_per_group )
{
	__local spinor scratch[128];
	__local double2 * d_scratch = (__local double2 *) scratch;

	__global const double2 * d_in = (__global double2 *) in;
	__global const double2 * d_out = (__global double2 *) out;

	size_t local_id = get_local_id(0);
	size_t id = get_group_id(0) * threads_per_group + local_id;

	for(size_t group_offset = get_group_id(0) * 128; group_offset < elems; group_offset += get_global_size(0)) {
		// make sure you don't run off the end
		size_t elem_count = 128 * sizeof(spinor) / sizeof(double2);
		event_t copy_in_event = async_work_group_copy(d_scratch, &d_in[group_offset], elem_count, 0);
		event_t copy_out_event = async_work_group_copy(&d_out[group_offset], d_scratch, elem_count, copy_in_event);
		wait_group_events(1, copy_out_event);
	}
}
