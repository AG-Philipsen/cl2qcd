(* ::Package:: *)

BeginPackage["Vectors`"]

squareNorm::usage = 
	"squareNorm[v] takes a generic vector as argument and compute its norm squared."

Begin["Private`"]

squareNorm[v_]:=
	Module[{squareNorm=Norm[v]^2},
	squareNorm
	]

End[]

EndPackage[]
