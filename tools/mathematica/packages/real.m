(* ::Package:: *)

BeginPackage["Real`"]

nonTrivialRealPar::usage = 
	"nonTrivialRealPar gives the real number 0.123456"

Begin["Private`"]

nonTrivialRealPar:=
	Module[ {nonTrivialRealPar=0.123456},
	nonTrivialRealPar
	]

End[]

EndPackage[]
