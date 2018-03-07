/**
 * Copyright (c) 2012 Matthias Bach
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

/*******************************************************************************
*
* file ranlxd.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

#ifndef RANLXD_H
#define RANLXD_H

extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

#endif
