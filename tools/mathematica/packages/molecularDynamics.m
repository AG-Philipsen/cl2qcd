(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["real.m",Path->{NotebookDirectory[]}]
Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["threeComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["eightComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["vectors.m",Path->{NotebookDirectory[]}]
Get["gammaMatricesTensorColorIdentity.m",Path->{NotebookDirectory[]}]


(*Forces for a specific site,for constant gaugemomenta, gauge- and \
spinorfield*)


BeginPackage["Forces`"]

ForceGauge::usage = 
 	""

ForceGaugeTlsym::usage = 
 	""

LinkByDiracTrace::usage = 
 	""

ForceFermionWilson::usage = 
 	""

Begin["Private`"]
Needs["GammaMatricesTensorColorIdentity`"]
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

ForceFermionWilson[gm_, k_, s_, u_, bc_] :=
 	Module[ {ForceFermionWilson = squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma4[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma1[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma2[cold3x3mat],bc]]]
                                  + squareNorm[gm + algebraElement[LinkByDiracTrace[k,s,u,Gamma3[cold3x3mat],bc]]]},
  	ForceFermionWilson
  	]

End[]

EndPackage[]






