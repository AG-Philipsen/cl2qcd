(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["real.m",Path->{NotebookDirectory[]}]
Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["threeComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["vectors.m",Path->{NotebookDirectory[]}]


(*Fermionmatrices for a specific site,for constant gauge- and \
spinorfield*)


BeginPackage["StaggeredDiracUtilities`"]

etaBC::usage = 
 	"Factor modifying the staggered phases according to the used Boundary Conditions"

etaCP::usage = 
 	"Factor multiplying lonks in temporal direction according to the used (Imaginary) Chemical Potential"

Begin["Private`"]

etaBC[th_, n_] :=
 	Module[ {etaBC = (0.5*Cos[th*Pi/n] + I*0.5*Sin[th*Pi/n]) },
  	etaBC
  	]

etaCP[cp_] :=
 	Module[ {etaCP = (Cos[cp] + I*Sin[cp]) },
  	etaCP
  	]

End[]

EndPackage[]


BeginPackage["MStaggered`"]

chiPosPhasePosMu::usage = 
 	"three-components vector obtained by the application of the link in +mu-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a positive
		staggered phase for that site and direction"

chiNegPhasePosMu::usage = 
 	"three-components vector obtained by the application of the link in +mu-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a negative
		staggered phase for that site and direction"

chiPosPhaseNegMu::usage = 
 	"three-components vector obtained by the application of the link in -mu-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a positive
		staggered phase for that site and direction"

chiNegPhaseNegMu::usage = 
 	"three-components vector obtained by the application of the link in -mu-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a negative
		staggered phase for that site and direction"

MStaggeredBC::usage =
	"M normal staggered fermionmatrix (see fermionmatrix_staggered_M.cl)"

Begin["Private`"]
Needs["StaggeredDiracUtilities`"]
Needs["Vectors`"]

chiPosPhasePosMu[u_, s_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = etaBC[th, n] * u.s },
  	chiPosPhasePosMu
  	]

chiNegPhasePosMu[u_, s_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = -etaBC[th, n] * u.s },
  	chiPosPhasePosMu
  	]

chiPosPhaseNegMu[u_, s_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = Conjugate[etaBC[th, n]]*ConjugateTranspose[u].s },
  	chiPosPhasePosMu
  	]

chiNegPhaseNegMu[u_, s_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = Conjugate[-etaBC[th, n]]*ConjugateTranspose[u].s },
  	chiPosPhasePosMu
  	]

MStaggeredBC[u_, s_, m_, tht_, ths_, nt_, ns_] :=
 	Module[ {MStaggeredBC = ns^3*nt*((1/4)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      chiPosPhasePosMu[u, s, tht, nt] - 
      chiPosPhaseNegMu[u, s, tht, nt]] + (1/4)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      chiNegPhasePosMu[u, s, tht, nt] - 
      chiNegPhaseNegMu[u, s, tht, nt]] + (1/8)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      2*(chiNegPhasePosMu[u, s, ths, ns] - 
         chiNegPhaseNegMu[u, s, ths, ns]) + 
      chiPosPhasePosMu[u, s, tht, nt] - 
      chiPosPhaseNegMu[u, s, tht, nt]] + (1/8)*
    squareNorm[
     m*s + 3*(chiPosPhasePosMu[u, s, ths, ns] - 
         chiPosPhaseNegMu[u, s, ths, ns]) + 
      chiPosPhasePosMu[u, s, tht, nt] - 
      chiPosPhaseNegMu[u, s, tht, nt]] + (1/8)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      2*(chiNegPhasePosMu[u, s, ths, ns] - 
         chiNegPhaseNegMu[u, s, ths, ns]) + 
      chiNegPhasePosMu[u, s, tht, nt] - 
      chiNegPhaseNegMu[u, s, tht, nt]] + (1/8)*
    squareNorm[
     m*s + 3*(chiPosPhasePosMu[u, s, ths, ns] - 
         chiPosPhaseNegMu[u, s, ths, ns]) + 
      chiNegPhasePosMu[u, s, tht, nt] - 
      chiNegPhaseNegMu[u, s, tht, nt]])},
  	MStaggeredBC
  	]

End[]

EndPackage[]





BeginPackage["DKS`"]

chiPosPhasePosT::usage = 
 	"three-components vector obtained by the application of the link in +t-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a positive
		staggered phase for that site and direction"

chiNegPhasePosT::usage =
 	"three-components vector obtained by the application of the link in +t-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a negative
		staggered phase for that site and direction"

chiPosPhaseNegT::usage = 
 	"three-components vector obtained by the application of the link in -t-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a positive
		staggered phase for that site and direction"

chiNegPhaseNegT::usage = 
 	"three-components vector obtained by the application of the link in -t-direction to the
		spinor on the local site and by multiplication with the BC factor assuming a negative
		staggered phase for that site and direction"

DKStaggered::usage =
	"Even-odd or Odd-even block of the staggered standard Dirac operator (see fermionmatrix_staggered_eo_DKS.cl)"

Begin["Private`"]
Needs["StaggeredDiracUtilities`"]
Needs["Vectors`"]
Needs["MStaggered`"]

chiPosPhasePosT[u_, s_, cp_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = etaBC[th, n]* etaCP[cp] * u.s },
  	chiPosPhasePosMu
  	]

chiNegPhasePosT[u_, s_, cp_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = -etaBC[th, n]* etaCP[cp] * u.s },
  	chiPosPhasePosMu
  	]

chiPosPhaseNegT[u_, s_, cp_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = Conjugate[etaBC[th, n]] * ConjugateTranspose[etaCP[cp] * u].s },
  	chiPosPhasePosMu
  	]

chiNegPhaseNegT[u_, s_, cp_, th_, n_] :=
 	Module[ {chiPosPhasePosMu = Conjugate[-etaBC[th, n]] * ConjugateTranspose[etaCP[cp] * u].s },
  	chiPosPhasePosMu
  	]

DKStaggered[u_, s_, m_, cp_, tht_, ths_, nt_, ns_]:=
	Module[ {DKStaggered = ns^3*nt*((1/4)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      chiPosPhasePosT[u, s, cp, tht, nt] - 
      chiPosPhaseNegT[u, s, cp, tht, nt]] + (1/4)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      chiNegPhasePosT[u, s, cp, tht, nt] - 
      chiNegPhaseNegT[u, s, cp, tht, nt]] + (1/8)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      2*(chiNegPhasePosMu[u, s, ths, ns] - 
         chiNegPhaseNegMu[u, s, ths, ns]) + 
      chiPosPhasePosT[u, s, cp, tht, nt] - 
      chiPosPhaseNegT[u, s, cp, tht, nt]] + (1/8)*
    squareNorm[
     m*s + 3*(chiPosPhasePosMu[u, s, ths, ns] - 
         chiPosPhaseNegMu[u, s, ths, ns]) + 
      chiPosPhasePosT[u, s, cp, tht, nt] - 
      chiPosPhaseNegT[u, s, cp, tht, nt]] + (1/8)*
    squareNorm[
     m*s + chiPosPhasePosMu[u, s, ths, ns] - 
      chiPosPhaseNegMu[u, s, ths, ns] + 
      2*(chiNegPhasePosMu[u, s, ths, ns] - 
         chiNegPhaseNegMu[u, s, ths, ns]) + 
      chiNegPhasePosT[u, s, cp, tht, nt] - 
      chiNegPhaseNegT[u, s, cp, tht, nt]] + (1/8)*
    squareNorm[
     m*s + 3*(chiPosPhasePosMu[u, s, ths, ns] - 
         chiPosPhaseNegMu[u, s, ths, ns]) + 
      chiNegPhasePosT[u, s, cp, tht, nt] - 
      chiNegPhaseNegT[u, s, cp, tht, nt]])},
	DKStaggered
	]

End[]

EndPackage[]




