(* ::Package:: *)

(*
 * Copyright 2012,2013 Lars Zeidlewicz,Christopher Pinke,
 * Matthias Bach,Christian Sch\[ADoubleDot]fer,Stefano Lottini,Alessandro Sciarra
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
*)


BeginPackage["GammaMatricesTensorColorIdentity`"]

Gamma1::usage =
	"Gamma1[x] gives the Gell-Mann matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Gamma2::usage =
	"Gamma2[x] gives the Gell-Mann matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Gamma3::usage =
	"Gamma3[x] gives the Gell-Mann matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Gamma4::usage =
	"Gamma3[x] gives the Gell-Mann matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Gamma5::usage =
	"Gamma5[x] gives the Gell-Mann matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

UnitMatrix::usage =
	"UnitMatrix[x] gives the Identity 4-by-4 matrix, tensor multiplied
	with the 3-by-3 x matrix in color space."

Begin["Private`"]

Gamma1[x_] :=
	Module[ {gamma1=ArrayFlatten[{{0, 0, 0, -I*x}, {0, 0, -I*x, 0}, {0, I*x, 0,
    0}, {I*x, 0, 0, 0}}]},
	gamma1
	]

Gamma2[x_] :=
	Module[ {gamma2=ArrayFlatten[{{0, 0, 0, -1.*x}, {0, 0, x*1., 0}, {0, 1.*x, 0,
   0}, {-1.*x, 0, 0, 0}}]},
	gamma2
	]

Gamma3[x_] :=
	Module[ {gamma3=ArrayFlatten[{{0, 0, -I*x, 0}, {0, 0, 0, x*I}, {I*x, 0, 0,
   0}, {0, -I*x, 0, 0}}]},
	gamma3
	]

Gamma4[x_] :=
	Module[ {gamma4=ArrayFlatten[{{0, 0, -1.*x, 0}, {0, 0, 0, x*-1.}, {-1.*x, 0, 0,
   0}, {0, -1.*x, 0, 0}}]},
	gamma4
	]

Gamma5[x_] :=
	Module[ {gamma5=ArrayFlatten[{{1.*x, 0, 0, 0}, {0, 1.*x, 0, 0}, {0, 0, -1.*x, 0}, {0, 0,
    0, -1.*x}}]},
	gamma5
	]

UnitMatrix[x_] :=
	Module[ {unitMatrix=ArrayFlatten[{{x, 0, 0, 0}, {0, x, 0, 0}, {0, 0, x, 0}, {0, 0, 0, x}}]},
	unitMatrix
	]

End[]

EndPackage[]
