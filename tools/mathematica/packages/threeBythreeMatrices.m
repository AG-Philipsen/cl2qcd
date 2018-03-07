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


BeginPackage["ThreeBythreeMatrices`"]

cold3x3mat::usage =
	"coldSU3mat gives the 3-by-3 Identity matrix"

nonTrivial3x3mat::usage =
	"nonTrivialSU3mat gives a non-trivially filled 3-by-3 matrix"

stapleSum::usage =
	"stapleSum[s] takes an SU(3) matrix as argument and compute the sum of staples."

rectStapleSum::usage =
	"rectStapleSum[s] takes an SU(3) matrix as argument and compute the sum of rectangles staples."

tracelessAntihermitianPart::usage =
	"tracelessAntihermitianPart takes a 3-by-3 matrix as argument and computes its traceless
	 anti-hermitian part which is another 3-by-3 matrix."

mat3x3FromKroneckerProductOf3ComponentsVectors::usage=
	"mat3x3FromKroneckerProductOf3ComponentsVectors calculates the Dirac-Trace of the matrix resulting
	from multiplying U*V^dagger =  u*v^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors
	(using spinprojection) see tr_v_times_u_dagger in operations_su3vec.cl"

Begin["Private`"]

cold3x3mat:=
	Module[ {coldSU3mat=IdentityMatrix[3]},
	coldSU3mat
	]

nonTrivial3x3mat:=
	Module[ {nonTrivialSU3mat={{0.130189 + I*0.260378, 0.260378 + I*0.390567,
  0.520756 + I*0.650945}, {.572742 + I*.403041, .371222 +
   I*.321726, -.449002 + I*(-.258088)}, {0.0000000000000000111022 +
   I*.651751, .0271563 + I*(-.733219), -.0271563 + I*.190094}}},
	nonTrivialSU3mat
	]

stapleSum[u_]:=
	Module[{stapleSum=3*(u.ConjugateTranspose[u].ConjugateTranspose[u] +
   ConjugateTranspose[u].ConjugateTranspose[u].u)},
	stapleSum
	]

rectStapleSum[u_]:=
	Module[{rectStapleSum=3*(u.u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u] +
							 u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u +
							 ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u.u +
							 ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u.u +
							 u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u].u +
							 u.u.ConjugateTranspose[u].ConjugateTranspose[u].ConjugateTranspose[u])},
	rectStapleSum
	]

tracelessAntihermitianPart[u_]:=
	Module[{tracelessAntihermitianPart= (1/2)*(u - ConjugateTranspose[u]) -
										(1/6)*Tr[u - ConjugateTranspose[u]]*IdentityMatrix[3]},
	tracelessAntihermitianPart
	]

mat3x3FromKroneckerProductOf3ComponentsVectors[s1_,s2_,s3_,s4_]:=
	Module[{mat3x3FromKroneckerProductOf3ComponentsVectors=KroneckerProduct[s1, ConjugateTranspose[s2]]
															+ KroneckerProduct[s3, ConjugateTranspose[s4]]},
	mat3x3FromKroneckerProductOf3ComponentsVectors
	]

End[]

EndPackage[]
