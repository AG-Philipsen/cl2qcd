(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["gammaMatricesTensorColorIdentity.m",Path->{NotebookDirectory[]}]
Get["real.m",Path->{NotebookDirectory[]}]
Get["threeBythreeMatrices.m",Path->{NotebookDirectory[]}]
Get["twelveComponentsVectors.m",Path->{NotebookDirectory[]}]
Get["vectors.m",Path->{NotebookDirectory[]}]


(*Fermionmatrices for a specific site,for constant gauge- and \
spinorfield*)


BeginPackage["DslashWilson`"]

Dslash::usage = 
	"compute the dslash taking kappa, a spinor filling and a 3x3matrix filling"

DslashWithSpecificBoundaryConditions::usage = 
	"compute the dslash taking kappa, a spinor filling, a 3x3matrix filling and
		taking into account specific boundary conditions both in temporal and spatial directions"

DslashWithSpecificBoundaryConditionsAndImaginaryChemicalPotential::usage = 
	"compute the dslash taking kappa, a spinor filling, a 3x3matrix filling,
		taking into account specific boundary conditions both in temporal and spatial directions
		and at a given imaginary value for the chemical potential"

DslashWithSpecificBoundaryConditionsAndChemicalPotential::usage = 
	"compute the dslash taking kappa, a spinor filling, a 3x3matrix filling,
		taking into account specific boundary conditions both in temporal and spatial directions
		and at a given complex chemical potential"

DslashWithImaginaryChemicalPotential::usage = 
	"compute the dslash taking kappa, a spinor filling, a 3x3matrix filling
		and at a given imaginary value for the chemical potential"

DslashWithRealChemicalPotential::usage = 
	"compute the dslash taking kappa, a spinor filling, a 3x3matrix filling
		and at a given real value for the chemical potential"

Begin["Private`"]
Needs["GammaMatricesTensorColorIdentity`"]

Dslash[k_, s_, u_]:=
	Module[ {Dslash=k*((UnitMatrix[u] - Gamma1[u]).s
					+ (UnitMatrix[u] - Gamma2[u]).s
					+ (UnitMatrix[u] - Gamma3[u]).s
					+ (UnitMatrix[u] - Gamma4[u]).s
					+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma1[u\[ConjugateTranspose]]).s
					+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma2[u\[ConjugateTranspose]]).s
					+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma3[u\[ConjugateTranspose]]).s
					+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma4[u\[ConjugateTranspose]]).s)},
	Dslash
	]

DslashWithSpecificBoundaryConditions[k_, s_, u_, ns_, nt_, ths_, tht_]:=
	Module[ {DslashWithSpecificBoundaryConditions=k*((UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma1[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma2[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma3[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[tht * I * Pi / nt]] - Gamma4[u * Exp[tht * I * Pi / nt]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma1[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma2[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma3[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]] + Gamma4[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]]).s)},
	DslashWithSpecificBoundaryConditions
	]

DslashWithSpecificBoundaryConditionsAndImaginaryChemicalPotential[k_, s_, u_, ns_, nt_, ths_, tht_, cp_]:=
	Module[ {DslashWithSpecificBoundaryConditionsAndImaginaryChemicalPotential=k*((UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma1[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma2[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma3[u * Exp[ths * I * Pi / ns]]).s
													+ Exp[+ I * cp] * (UnitMatrix[u * Exp[tht * I * Pi / nt]] - Gamma4[u * Exp[tht * I * Pi / nt]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma1[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma2[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma3[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ Exp[+ I * cp] * (UnitMatrix[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]] + Gamma4[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]]).s)},
	DslashWithSpecificBoundaryConditionsAndImaginaryChemicalPotential
	]

DslashWithSpecificBoundaryConditionsAndChemicalPotential[k_, s_, u_, ns_, nt_, ths_, tht_, cpRe_, cpIm_]:=
	Module[ {DslashWithSpecificBoundaryConditionsAndChemicalPotential=k*((UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma1[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma2[u * Exp[ths * I * Pi / ns]]).s
													+ (UnitMatrix[u * Exp[ths * I * Pi / ns]] - Gamma3[u * Exp[ths * I * Pi / ns]]).s
													+ Exp[+ I * cpIm] * Exp[+ cpRe] * (UnitMatrix[u * Exp[tht * I * Pi / nt]] - Gamma4[u * Exp[tht * I * Pi / nt]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma1[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma2[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ (UnitMatrix[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]] + Gamma3[u\[ConjugateTranspose] * Exp[- ths * I * Pi / ns]]).s
													+ Exp[- I * cpIm] * Exp[- cpRe] * (UnitMatrix[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]] + Gamma4[u\[ConjugateTranspose] * Exp[- tht * I * Pi / nt]]).s)},
	DslashWithSpecificBoundaryConditionsAndChemicalPotential
	]

DslashWithImaginaryChemicalPotential[k_, s_, u_, cpIm_]:=
	Module[ {DslashWithImaginaryChemicalPotential=k * ((UnitMatrix[u] - Gamma1[u]).s
													+ (UnitMatrix[u] - Gamma2[u]).s
													+ (UnitMatrix[u] - Gamma3[u]).s
													+ Exp[+ I * cpIm] * (UnitMatrix[u] - Gamma4[u]).s
													+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma1[u\[ConjugateTranspose]]).s
													+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma2[u\[ConjugateTranspose]]).s
													+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma3[u\[ConjugateTranspose]]).s
													+ Exp[- I * cpIm] * (UnitMatrix[u\[ConjugateTranspose]] + Gamma4[u\[ConjugateTranspose]]).s)},
	DslashWithImaginaryChemicalPotential
	]

DslashWithRealChemicalPotential[k_, s_, u_, cpRe_]:=
	Module[ {DslashWithRealChemicalPotential=k * ((UnitMatrix[u] - Gamma1[u]).s
												+ (UnitMatrix[u] - Gamma2[u]).s
												+ (UnitMatrix[u] - Gamma3[u]).s
												+ Exp[cpRe] * (UnitMatrix[u] - Gamma4[u]).s
												+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma1[u\[ConjugateTranspose]]).s
												+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma2[u\[ConjugateTranspose]]).s
												+ (UnitMatrix[u\[ConjugateTranspose]] + Gamma3[u\[ConjugateTranspose]]).s
												+ Exp[cpRe] * (UnitMatrix[u\[ConjugateTranspose]] + Gamma4[u\[ConjugateTranspose]]).s)},
	DslashWithRealChemicalPotential
	]

End[]

EndPackage[]


BeginPackage["MWilson`"]

MWilson::usage = 
	"compute the Wilson matrix taking kappa, a spinor filling and a 3x3matrix filling"

Begin["Private`"]
Needs["DslashWilson`"]

MWilson[k_, s_, u_]:=
	Module[ {MWilson = s - Dslash[k, s, u]},
	MWilson
	]

End[]

EndPackage[]


BeginPackage["MTwistedMass`"]

MTwistedMassPlus::usage = 
	""

MTwistedMassMinus::usage = 
	""

MTwistedMassSitediagonal::usage = 
	""

MTwistedMassSitediagonalMinus::usage = 
	""

MTwistedMassInverseSitediagonal::usage = 
	""

MTwistedMassInverseSitediagonalMinus::usage = 
	""

MTwistedMassSitediagonalAndGamma5EvenOdd::usage = 
	""

MTwistedMassSitediagonalMinusAndGamma5EvenOdd::usage = 
	""

Begin["Private`"]
Needs["DslashWilson`"]
Needs["GammaMatricesTensorColorIdentity`"]

MTwistedMassPlus[k_, s_, u_, m_]:=
	Module[ {MTwistedMassPlus = (IdentityMatrix[12] + 2. I*k*m*Gamma5[IdentityMatrix[3]]).s - 
 Dslash[k, s, u]},
	MTwistedMassPlus
	]

MTwistedMassMinus[k_, s_, u_, m_]:=
	Module[ {MTwistedMassMinus = (IdentityMatrix[12] - 2.*I*k*m*Gamma5[IdentityMatrix[3]]).s - 
 Dslash[k, s, u]},
	MTwistedMassMinus
	]

MTwistedMassSitediagonal[s_, k_, m_]:=
	Module[ {MTwistedMassSitediagonal = (IdentityMatrix[12] + 2.*I*m*k*Gamma5[IdentityMatrix[3]]).s},
	MTwistedMassSitediagonal
	]

MTwistedMassSitediagonalMinus[s_, k_, m_]:=
	Module[ {MTwistedMassSitediagonalMinus = (IdentityMatrix[12] - 2.*I*m*k*Gamma5[IdentityMatrix[3]]).s},
	MTwistedMassSitediagonalMinus
	]

MTwistedMassInverseSitediagonal[s_, k_, m_]:=
	Module[ {MTwistedMassInverseSitediagonal = (1./(1. + 4*k*k*m*m))*(IdentityMatrix[12] - 
    2.*I*m*k*Gamma5[IdentityMatrix[3]]).s},
	MTwistedMassInverseSitediagonal
	]

MTwistedMassInverseSitediagonalMinus[s_, k_, m_]:=
	Module[ {MTwistedMassInverseSitediagonalMinus = (1./(1. + 4*k*k*m*m))*(IdentityMatrix[12] + 
    2.*I*m*k*Gamma5[IdentityMatrix[3]]).s},
	MTwistedMassInverseSitediagonalMinus
	]

MTwistedMassSitediagonalAndGamma5EvenOdd[s_, k_, m_]:=
	Module[ {MTwistedMassSitediagonalAndGamma5EvenOdd = Gamma5[IdentityMatrix[3]].MTwistedMassSitediagonal[s, k, m]},
	MTwistedMassSitediagonalAndGamma5EvenOdd
	]

MTwistedMassSitediagonalMinusAndGamma5EvenOdd[s_, k_, m_]:=
	Module[ {MTwistedMassSitediagonalMinusAndGamma5EvenOdd = Gamma5[IdentityMatrix[3]].MTwistedMassSitediagonalMinus[s, k, m]},
	MTwistedMassSitediagonalMinusAndGamma5EvenOdd
	]

End[]

EndPackage[]


BeginPackage["DslashTwistedMass`"]

DslashEvenOddAndMTmInverseSitediagonalMinus::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonal::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonalBC::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonalChemPotIm::usage = 
	""

DslashEvenOddAndMTmInverseSitediagonalChemPotRe::usage = 
	""

Begin["Private`"]
Needs["MTwistedMass`"]
Needs["DslashWilson`"]

DslashEvenOddAndMTmInverseSitediagonalMinus[s_, u_, k_, m_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalMinus = MTwistedMassInverseSitediagonalMinus[(-1.)*Dslash[k, s, u], k, m]},
	DslashEvenOddAndMTmInverseSitediagonalMinus
	]

DslashEvenOddAndMTmInverseSitediagonal[s_, u_, k_, m_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonal = MTwistedMassInverseSitediagonal[(-1.)*Dslash[k, s, u], k, m]},
	DslashEvenOddAndMTmInverseSitediagonal
	]


DslashEvenOddAndMTmInverseSitediagonalBC[s_, u_, k_, m_, ns_, nt_, ths_, tht_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalBC = MTwistedMassInverseSitediagonal[(-1.)*
  DslashWithSpecificBoundaryConditions[k, s, u, ns, nt, ths, tht], k, m]},
	DslashEvenOddAndMTmInverseSitediagonalBC
	]


DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[s_, u_, k_, m_, ns_, nt_, ths_, tht_, cp_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot = MTwistedMassInverseSitediagonal[(-1.)*
  DslashWithSpecificBoundaryConditionsAndImaginaryChemicalPotential[k, s, u, ns, nt, ths, tht, cp], k, m]},
	DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot
	]


DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[s_, u_, k_,m_, ns_, nt_, ths_, tht_, cpIm_, cpRe_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe = MTwistedMassInverseSitediagonal[(-1.)*
  DslashWithSpecificBoundaryConditionsAndChemicalPotential[k, s, u, ns, nt, ths, tht, cpIm, cpRe], k,
  m]},
	DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe
	]


DslashEvenOddAndMTmInverseSitediagonalChemPotIm[s_, u_, k_, m_, ns_,nt_, cpIm_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalChemPotIm = MTwistedMassInverseSitediagonal[(-1.)*
  DslashWithImaginaryChemicalPotential[k, s, u, cpIm], k, m]},
	DslashEvenOddAndMTmInverseSitediagonalChemPotIm
	]


DslashEvenOddAndMTmInverseSitediagonalChemPotRe[s_, u_, k_, m_, ns_, nt_, cpRe_]:=
	Module[ {DslashEvenOddAndMTmInverseSitediagonalChemPotRe = MTwistedMassInverseSitediagonal[(-1.)*
  DslashWithRealChemicalPotential[k, s, u, cpRe], k, m]},
	DslashEvenOddAndMTmInverseSitediagonalChemPotRe
	]


End[]

EndPackage[]
