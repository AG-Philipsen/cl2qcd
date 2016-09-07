(* ::Package:: *)

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
