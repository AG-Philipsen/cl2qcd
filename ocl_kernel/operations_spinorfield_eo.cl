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

/** @file
 * Device code for operations on the even-odd preconditioned spinor field
 */

// TODO document
inline spinor getSpinor_eo(__global const spinorStorageType* const restrict in, const uint idx)
{
#ifdef _USE_SOA_
    return (spinor){{// su3vec = 3 * cplx
                     in[0 * EOPREC_SPINORFIELD_STRIDE + idx], in[1 * EOPREC_SPINORFIELD_STRIDE + idx],
                     in[2 * EOPREC_SPINORFIELD_STRIDE + idx]},
                    {// su3vec = 3 * cplx
                     in[3 * EOPREC_SPINORFIELD_STRIDE + idx], in[4 * EOPREC_SPINORFIELD_STRIDE + idx],
                     in[5 * EOPREC_SPINORFIELD_STRIDE + idx]},
                    {// su3vec = 3 * cplx
                     in[6 * EOPREC_SPINORFIELD_STRIDE + idx], in[7 * EOPREC_SPINORFIELD_STRIDE + idx],
                     in[8 * EOPREC_SPINORFIELD_STRIDE + idx]},
                    {// su3vec = 3 * cplx
                     in[9 * EOPREC_SPINORFIELD_STRIDE + idx], in[10 * EOPREC_SPINORFIELD_STRIDE + idx],
                     in[11 * EOPREC_SPINORFIELD_STRIDE + idx]}};
#else
    return in[idx];
#endif
}

// TODO document
inline void putSpinor_eo(__global spinorStorageType* const restrict out, const uint idx, const spinor val)
{
#ifdef _USE_SOA_
    // su3vec = 3 * cplx
    out[0 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e0;
    out[1 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e1;
    out[2 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e2;

    // su3vec = 3 * cplx
    out[3 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e0;
    out[4 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e1;
    out[5 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e2;

    // su3vec = 3 * cplx
    out[6 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e0;
    out[7 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e1;
    out[8 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e2;

    // su3vec = 3 * cplx
    out[9 * EOPREC_SPINORFIELD_STRIDE + idx]  = val.e3.e0;
    out[10 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e1;
    out[11 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e2;
#else
    out[idx] = val;
#endif
}
