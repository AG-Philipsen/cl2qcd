(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/molecularDynamics.m"}], Path -> {NotebookDirectory[]}]


(* FGauge *)
ForceGauge[5.69,cold3x3mat,gaugeMomOne]
ForceGauge[5.69,nonTrivial3x3mat,gaugeMomOne]
ForceGauge[5.69,cold3x3mat,gaugeMomAscending]
ForceGauge[5.69,nonTrivial3x3mat,gaugeMomAscending]


(* FGaugeTlsym *)
c1:=-0.083333333
ForceGaugeTlsym[5.69,c1,cold3x3mat,gaugeMomOne]
ForceGaugeTlsym[5.69,c1,nonTrivial3x3mat,gaugeMomOne]
ForceGaugeTlsym[5.69,c1,cold3x3mat,gaugeMomAscending]
ForceGaugeTlsym[5.69,c1,nonTrivial3x3mat,gaugeMomAscending]


ForceFermionWilson[gaugeMomOne,nonTrivialRealPar,spinorCold,IdentityMatrix[3],(1+0*I)]
FermionForce[gaugeMomOne,0.123456,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]


FermionForce[gaugeMomAsc,0.123456,spinorCold,IdentityMatrix[3],(1+0*I)]
FermionForce[gaugeMomAsc,0.123456,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]


FermionForce[gaugeMomOne,0.123456,spinorCold,nonTrivialMatrix,(1+0*I)]
FermionForce[gaugeMomOne, 0.123456,spinorAscendingComplex,nonTrivialMatrix,(1+0*I)]


FermionForce[gaugeMomAsc,0.123456,spinorCold,nonTrivialMatrix,(1+0*I)]
FermionForce[gaugeMomAsc,0.123456,spinorAscendingComplex,nonTrivialMatrix,(1+0*I)]



