(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["real.m",Path->{NotebookDirectory[]}]
Get["gammaMatricesTensorColorIdentity.m",Path->{NotebookDirectory[]}]
Get["vectors.m",Path->{NotebookDirectory[]}]
Get["threeComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["eightComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]


(*Forces for a specific site,for constant gaugemomenta, gauge- and spinorfield*)


BeginPackage["Forces`"]

ForceGauge::usage = 
 	"ForceGauge performs the standard Wilson-action gauge force calculation given a gaugefield filling, a gaugemomentum filling and a value for the coupling beta"

ForceGaugeTlsym::usage = 
 	"ForceGaugeTlsym performs the tlSym improved gauge action gauge force calculation given a gaugefield filling, a gaugemomentum filling, a value for the coupling
      beta, and a value for the c1 coefficient in tlSym gauge action"

LinkByDiracTrace::usage = 
 	"LinkByDiracTrace performs the product, for a given direction \mu, of the local link in that direction by the trace in Dirac space of the cartesian product in
      color space of the appropriate spinors in the Wilson hopping term (takes into account both +mu and -mu).
      g is the generic gamma matrix tensor product with the col.identity"

LinkByDiracTraceEvenOrOdd::usage = 
 	"LinkByDiracTraceEvenOrOdd performs the same operation as LinkByDiracTrace, but only on the even/odd sublattice. Relevant for when the fermion_force_eo kernel is called with a boolean parameter EVEN/ODD.
      As a result only even/odd sites are considered and the force receives contribution (the same) only by either positive or negative mu direction computation. 
      The kernel must be called twice with different evenodd value to obtain the force on all sites and there is a dedicated test to check that this gives the full force of the nonEvenOdd case.
	  g is the generic gamma matrix tensor product with the col.identity"

ForceFermionWilson::usage = 
 	"ForceFermionWilson computes the local fermionic force for the Wilson case given the gaugemomenta field gm, the hopping parameter k, the spinorfield s, the gaugefield u, and the boundary conditions bc.
      It then performs direction-by-direction the squarenorm of the updated gaugemomenta and sum up over the 4 directions."

ForceFermionWilsonEvenOrOdd::usage = 
 	"ForceFermionWilson computes the local fermionic force for the Wilson case given the gaugemomenta field gm, the hopping parameter k, the spinorfield s, the gaugefield u, and the boundary conditions bc.
      It then performs direction-by-direction the squarenorm of the updated gaugemomenta and sum up over the 4 directions. Relevant for when the fermion_force_eo kernel is called with a boolean parameter EVEN/ODD.
      As a result only even/odd sites are considered and the force receives contribution (the same) only by either positive or negative mu direction computation."

ForceFermionStaggeredEven::usage = 
 	"ForceFermionStaggeredEven computes the local fermionic force for the Staggered case on even sites, given the gaugemomenta field gm, the gaugefield u, the staggered phase eta and the spinorfield s."

ForceFermionStaggeredOdd::usage = 
 	"ForceFermionStaggeredOdd computes the local fermionic force for the Staggered case on odd sites, given the gaugemomenta field gm, the gaugefield u, the staggered phase eta and the spinorfield s."

Begin["Private`"]
Needs["GammaMatricesTensorColorIdentity`"]
Needs["GellMannMatrices`"]
Needs["TwelveComponentsVectors`"]
Needs["EightComponentsVectors`"]
Needs["ThreeComponentsVectors`"]
Needs["ThreeBythreeMatrices`"]
Needs["Vectors`"]

ss[s_,g_]:=idtyPlusGenericGammaTimesSpinor[s,g]
g5ss[s_,g_]:=idtyPlusGenericGammaTimesSpinor[gamma5TimesSpinor[s],g]

ForceGauge[b_,u_,gm_] :=
 	Module[ {ForceGauge = squareNorm[- b / 3 * algebraElement[u.stapleSum[u]] + gm]},
  	ForceGauge
  	]

ForceGaugeTlsym[b_,c1_,u_,gm_] :=
 	Module[ {ForceGaugeTlsym = squareNorm[-c1 * b / 3 * algebraElement[u.rectStapleSum[u]] + gm]},
  	ForceGaugeTlsym
  	]

LinkByDiracTrace[k_,s_,u_,g_,bc_]:=
	Module[ {LinkByDiracTrace = 2*k*bc*(u.ConjugateTranspose[
								mat3x3FromKroneckerProductOf3ComponentsVectors[extractFirstThreeComponentsOfWilsonSpinor[ss[s,g]],extractFirstThreeComponentsOfWilsonSpinor[g5ss[s,g]],extractSecondThreeComponentsOfWilsonSpinor[ss[s,g]],extractSecondThreeComponentsOfWilsonSpinor[g5ss[s,g]]]]
								+u.ConjugateTranspose[
								mat3x3FromKroneckerProductOf3ComponentsVectors[extractFirstThreeComponentsOfWilsonSpinor[g5ss[s,-g]],extractFirstThreeComponentsOfWilsonSpinor[ss[s,-g]],extractSecondThreeComponentsOfWilsonSpinor[g5ss[s,-g]],extractSecondThreeComponentsOfWilsonSpinor[ss[s,-g]]]])},
	LinkByDiracTrace
	]

LinkByDiracTraceEvenOrOdd[k_,s_,u_,g_,bc_]:=
	Module[ {LinkByDiracTraceEvenOrOdd = 2*k*bc*(u.ConjugateTranspose[
								mat3x3FromKroneckerProductOf3ComponentsVectors[extractFirstThreeComponentsOfWilsonSpinor[g5ss[s,-g]],extractFirstThreeComponentsOfWilsonSpinor[ss[s,-g]],extractSecondThreeComponentsOfWilsonSpinor[g5ss[s,-g]],extractSecondThreeComponentsOfWilsonSpinor[ss[s,-g]]]])},
	LinkByDiracTraceEvenOrOdd
	]

ForceFermionWilson[gm_, k_, s_, u_, bc_] :=
 	Module[ {ForceFermionWilson = squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma4[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma1[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma2[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma3[cold3x3mat],bc]]]},
  	ForceFermionWilson
  	]

ForceFermionWilsonEvenOrOdd[gm_, k_, s_, u_, bc_] :=
 	Module[ {ForceFermionWilsonEvenOrOdd = squareNorm[gm + algebraElement[LinkByDiracTraceEvenOrOdd[k,s,u,Gamma4[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTraceEvenOrOdd[k,s,u,Gamma1[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTraceEvenOrOdd[k,s,u,Gamma2[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTraceEvenOrOdd[k,s,u,Gamma3[cold3x3mat],bc]]]},
  	ForceFermionWilsonEvenOrOdd
  	]

ForceFermionStaggeredEven[gm_, u_, eta_, s_] :=
 	Module[ {ForceFermionStaggeredEven = squareNorm[gm + aeFromSu3[-I*tracelessAntihermitianPart[u.KroneckerProduct[(eta*s),ConjugateTranspose[s]]]]]},
  	ForceFermionStaggeredEven
  	]

ForceFermionStaggeredOdd[gm_, u_, eta_, s_] :=
 	Module[ {ForceFermionStaggeredOdd = squareNorm[gm - aeFromSu3[-I*tracelessAntihermitianPart[u.KroneckerProduct[(eta*s),ConjugateTranspose[s]]]]]},
  	ForceFermionStaggeredOdd
  	]

End[]

EndPackage[]



