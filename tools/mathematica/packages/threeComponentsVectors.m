(* ::Package:: *)

BeginPackage["ThreeComponentsVectors`"]

su3vecZero::usage =
	"su3vecZero gives a zero 3-components vector (uniformly filled by {0.,0.})."

su3vecCold::usage =
	"su3vecCold gives a cold 3-components vector (uniformly filled by {1.,0.})."

su3vecAscendingReal::usage =
	"su3vecAscendingReal gives a 3-components vector filled with ascending real numbers."

su3vecAscendingComplex::usage =
	"su3vecAscendingComplex gives a zero 3-components vector filled with ascending complex numbers."

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

countSSf[s_]:=
	Module[{countSf= Re[s[[1]]] + Im[s[[1]]] + Re[s[[2]]] + Im[s[[2]]] + Re[s[[3]]] + Im[s[[3]]]},
	countSf
	]

End[]

EndPackage[]
