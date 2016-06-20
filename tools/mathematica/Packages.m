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

BeginPackage["GellMannMatrices`"]

Lambda1::usage = 
	"Lambda1 gives the 1st Gell-Mann matrix."

Lambda2::usage = 
	"Lambda2 gives the 2nd Gell-Mann matrix."

Lambda3::usage = 
	"Lambda3 gives the 3rd Gell-Mann matrix."

Lambda4::usage = 
	"Lambda4 gives the 4th Gell-Mann matrix."

Lambda5::usage = 
	"Lambda5 gives the 5th Gell-Mann matrix."

Lambda6::usage = 
	"Lambda6 gives the 6th Gell-Mann matrix."

Lambda7::usage = 
	"Lambda7 gives the 7th Gell-Mann matrix."

Lambda8::usage = 
	"Lambda8 gives the 8th Gell-Mann matrix."

Begin["Private`"]

Lambda1:=
	Module[ {Lambda1={{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}},
	Lambda1
	]

Lambda2:=
	Module[ {Lambda2={{0, -I, 0}, {I, 0, 0}, {0, 0, 0}}},
	Lambda2
	]

Lambda3:=
	Module[ {Lambda3={{1, 0, 0}, {0, -1, 0}, {0, 0, 0}}},
	Lambda3
	]

Lambda4:=
	Module[ {Lambda4={{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}},
	Lambda4
	]

Lambda5:=
	Module[ {Lambda5={{0, 0, -I}, {0, 0, 0}, {I, 0, 0}}},
	Lambda5
	]

Lambda6:=
	Module[ {Lambda6={{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}},
	Lambda6
	]

Lambda7:=
	Module[ {Lambda7={{0, 0, 0}, {0, 0, -I}, {0, I, 0}}},
	Lambda7
	]

Lambda8:=
	Module[ {Lambda8={{1, 0, 0}, {0, 1, 0}, {0, 0, -2}}/Sqrt[3]},
	Lambda8
	]

End[]

EndPackage[]

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

BeginPackage["TwelveComponentsVectors`"]

spinorCold::usage = 
	"spinorCold gives a cold filled spinor (uniformly filled by {1.,0.})."

spinorAscendingReal::usage = 
	"spinorAscendingReal gives a spinor filled with ascending real numbers."

spinorAscendingComplex::usage = 
	"spinorAscendingReal gives a spinor filled with ascending complex numbers."

countSf::usage = 
	"countSf[s] takes a spinor as argument and compute the sum of all its (real+imaginary) components."

gamma5TimesSpinor::usage = 
	"gamma5TimesSpinor takes a spinor and applies Gamma5[Idty] to it."

Begin["Private`"]
Needs["GammaMatricesTensorColorIdentity`"]
f[x_, y_] := x; 
g[x_, y_] := (2*(x - 1) + 1)*1.0000000000000 + I*(2*(x - 1) + 2);

spinorCold := 
	Module[ {spinorCold=ConstantArray[1., {12, 1}]},
	spinorCold
	]

spinorAscendingReal:=
	Module[ {spinorAscendingReal=Array[ f, {12, 1}]},
	spinorAscendingReal
	]

spinorAscendingComplex:=
	Module[ {spinorAscendingComplex=Array[ g, {12, 1}]},
	spinorAscendingComplex
	]

countSf[s_]:=
	Module[{countSf = Re[s[[1]]] + Im[s[[1]]] + Re[s[[2]]] + Im[s[[2]]] + Re[s[[3]]] + Im[s[[3]]]  
                     + Re[s[[4]]] + Im[s[[4]]] + Re[s[[5]]] + Im[s[[5]]] + Re[s[[6]]] + Im[s[[6]]] 
                     + Re[s[[7]]] + Im[s[[7]]] + Re[s[[8]]] + Im[s[[8]]] + Re[s[[9]]] + Im[s[[9]]] 
                     + Re[s[[10]]] + Im[s[[10]]] + Re[s[[11]]] + Im[s[[11]]] + Re[s[[12]]] + Im[s[[12]]]},
	countSf
	]

gamma5TimesSpinor[s_]:=
	Module[{gamma5TimesSpinor=Gamma5[IdentityMatrix[3]].s},
	gamma5TimesSpinor
	]

End[]

EndPackage[]

BeginPackage["EightComponentsVectors`"]

gaugeMomZero::usage = 
	"gaugeMomZero gives a zero algebra vector (uniformly filled by {0.})."

gaugeMomOne::usage = 
	"gaugeMomOne gives a cold algebra vector (uniformly filled by {1.})."

gaugeMomAscending::usage =
	"gaugeMomAscending gives an algebra vector filled with ascending numbers."

algebraElement::usage =
	"algebraElement takes an SU(3) matrix as argument and computes the algebra element
	 Tr[i * T_k * (in - in^\\dag)], i.e. calculates the trace of I times a generator of the
	 algebra su(3) times a 3x3-matrix."

aeFromSu3::usage = 
	"aeFromSu3 takes an SU(3) matrix and builds an algebra element out of it."

Begin["Private`"]
Needs["GellMannMatrices`"]

gaugeMomZero:=
	Module[ {gaugeMomZero={0., 0., 0., 0., 0., 0., 0., 0.}},
	gaugeMomZero
	]

gaugeMomOne:=
	Module[ {gaugeMomOne={1., 1., 1., 1., 1., 1., 1., 1.}},
	gaugeMomOne
	]

gaugeMomAscending:=
	Module[ {gaugeMomAsc={1., 2., 3., 4., 5., 6., 7., 8.}},
	gaugeMomAsc
	]

algebraElement[u_]:=
	Module[{algebraElement={Tr[I*Lambda1.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda2.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda3.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda4.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda5.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda6.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda7.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda8.(u - ConjugateTranspose[u])]/2}},
	algebraElement
	]

aeFromSu3[u_]:=
	Module[{aeFromSu3={2*Re[u[[1, 2]]], -2*Im[u[[1, 2]]], Re[u[[1, 1]]] - Re[u[[2, 2]]],
						2*Re[u[[1, 3]]], -2 Im[u[[1, 3]]], 2*Re[u[[2, 3]]], -2*Im[u[[2, 3]]],
						(Re[u[[1, 1]]] + Re[u[[2, 2]]])*3/Sqrt[3]}},
	aeFromSu3
	]

End[]

EndPackage[]

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

End[]

EndPackage[]

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
