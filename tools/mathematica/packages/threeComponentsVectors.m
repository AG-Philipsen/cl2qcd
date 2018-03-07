(* ::Package:: *)

(*
 * Copyright 2012,2013 Lars Zeidlewicz,Christopher Pinke,
 * Matthias Bach,Christian Sch\[ADoubleDot]fer,Stefano Lottini,Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software:you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation,either version 3 of the License,or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY;without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.If not,see<http://www.gnu.org/licenses/>.
*)


BeginPackage["ThreeComponentsVectors`"]

su3vecZero::usage =
	"su3vecZero gives a zero 3-components vector (uniformly filled by {0.,0.})."

su3vecCold::usage =
	"su3vecCold gives a cold 3-components vector (uniformly filled by {1.,0.})."

su3vecAscendingReal::usage =
	"su3vecAscendingReal gives a 3-components vector filled with ascending real numbers."

su3vecAscendingComplex::usage =
	"su3vecAscendingComplex gives a zero 3-components vector filled with ascending complex numbers."

extractFirstThreeComponentsOfWilsonSpinor::usage =
	"extractFirstThreeComponentsOfWilsonSpinor build a 3-components vector using the first 3 entries of a 12 components one."

extractSecondThreeComponentsOfWilsonSpinor::usage =
	"extractFirstThreeComponentsOfWilsonSpinor build a 3-components vector using the second 3 entries of a 12 components one."

countSSf::usage =
	"countSSf[s] takes an su3vec as argument and compute the sum of all its (real+imaginary) components."

Begin["Private`"]
f[x_, y_] := x;
g[x_, y_] := (2*(x - 1) + 1)*1.0000000000000 + I*(2*(x - 1) + 2);

su3vecZero:=
	Module[ {su3vecZero=ConstantArray[0., {3, 1}]},
	su3vecZero
	]

su3vecCold:=
	Module[ {su3vecCold=ConstantArray[1., {3, 1}]},
	su3vecCold
	]

su3vecAscendingReal:=
	Module[ {su3vecAscendingReal=Array[f, {3, 1}]},
	su3vecAscendingReal
	]

su3vecAscendingComplex:=
	Module[ {su3vecAscendingComplex=Array[g, {3, 1}]},
	su3vecAscendingComplex
	]

extractFirstThreeComponentsOfWilsonSpinor[spinor_]:=
	Module[{extractFirstThreeComponentsOfWilsonSpinor={spinor[[1]],spinor[[2]],spinor[[3]]}},
	extractFirstThreeComponentsOfWilsonSpinor
	]

extractSecondThreeComponentsOfWilsonSpinor[spinor_]:=
	Module[{extractSecondThreeComponentsOfWilsonSpinor={spinor[[4]],spinor[[5]],spinor[[6]]}},
	extractSecondThreeComponentsOfWilsonSpinor
	]

countSSf[s_]:=
	Module[{countSf= Re[s[[1]]] + Im[s[[1]]] + Re[s[[2]]] + Im[s[[2]]] + Re[s[[3]]] + Im[s[[3]]]},
	countSf
	]

End[]

EndPackage[]
