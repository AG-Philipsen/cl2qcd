(* ::Package:: *)

Get["gammaMatricesTensorColorIdentity.m", Path -> {NotebookDirectory[]}]
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

idtyPlusGenericGammaTimesSpinor::usage = 
	"idtyPlusGenericGammaTimesSpinor takes a spinor and performs the operation (1 pm gamma_mu).spinor."

SaxpyAndGamma5EvenOdd::usage = 
	"SaxpyAndGamma5EvenOdd takes two spinors applies Gamma5[Idty] to the first and adds the result to the second."

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

idtyPlusGenericGammaTimesSpinor[s_,g_]:=
	Module[{idtyPlusGenericGammaTimesSpinor=(UnitMatrix[IdentityMatrix[3]]+g).s},
	idtyPlusGenericGammaTimesSpinor
	]

SaxpyAndGamma5EvenOdd[s1_, s2_, a_]:=
	Module[{SaxpyAndGamma5EvenOdd=Gamma5[IdentityMatrix[3]].(s2 - a*s1)},
	SaxpyAndGamma5EvenOdd
	]

End[]

EndPackage[]
