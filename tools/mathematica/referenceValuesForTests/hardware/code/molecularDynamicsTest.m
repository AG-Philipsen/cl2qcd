(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/molecularDynamics.m"}], Path -> {NotebookDirectory[]}]


FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/molecularDynamics.m"}]


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


(*force:Fermion Force*)
ForceFermionWilson[gaugeMomOne,nonTrivialRealPar,spinorCold,IdentityMatrix[3],(1+0*I)]
ForceFermionWilson[gaugeMomOne,nonTrivialRealPar,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]
ForceFermionWilson[gaugeMomAscending,nonTrivialRealPar,spinorCold,IdentityMatrix[3],(1+0*I)]
ForceFermionWilson[gaugeMomAscending,nonTrivialRealPar,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]
ForceFermionWilson[gaugeMomOne,nonTrivialRealPar,spinorCold,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilson[gaugeMomOne,nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilson[gaugeMomAscending,nonTrivialRealPar,spinorCold,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilson[gaugeMomAscending,nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,(1+0*I)]


(*force:Fermion Force EvenOdd.
In this case the kernel is called with a boolean parameter EVEN/ODD.
As a result only even/odd sites are considered and the force receives contribution (the same)
only by either positive or negative mu direction computation.
The kernel must be called twice with different evenodd value to obtain the force on all sites
and there is a dedicated test to check that this gives the full force of the nonEvenOdd case.*)
ForceFermionWilsonEvenOrOdd[gaugeMomOne,nonTrivialRealPar,spinorCold,IdentityMatrix[3],(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomOne,nonTrivialRealPar,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomAscending,nonTrivialRealPar,spinorCold,IdentityMatrix[3],(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomAscending,nonTrivialRealPar,spinorAscendingComplex,IdentityMatrix[3],(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomOne,nonTrivialRealPar,spinorCold,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomOne,nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomAscending,nonTrivialRealPar,spinorCold,nonTrivial3x3mat,(1+0*I)]
ForceFermionWilsonEvenOrOdd[gaugeMomAscending,nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,(1+0*I)]


(*force:Fermion Force Staggered with gaugeMomOne and nonTrivial3x3mat*)
(*Even, su3vecCold*)
ForceFermionStaggeredEven[gaugeMomOne,nonTrivial3x3mat,1.,su3vecCold]
ForceFermionStaggeredEven[gaugeMomOne,nonTrivial3x3mat,-1.,su3vecCold]
(*Odd, su3vecCold*)
ForceFermionStaggeredOdd[gaugeMomOne,nonTrivial3x3mat,1.,su3vecCold]
ForceFermionStaggeredOdd[gaugeMomOne,nonTrivial3x3mat,-1.,su3vecCold]
(*Even, su3vecAscendingComplex*)
ForceFermionStaggeredEven[gaugeMomOne,nonTrivial3x3mat,+1.,su3vecAscendingComplex]
ForceFermionStaggeredEven[gaugeMomOne,nonTrivial3x3mat,-1.,su3vecAscendingComplex]
(*Odd, su3vecAscendingComplex*)
ForceFermionStaggeredOdd[gaugeMomOne,nonTrivial3x3mat,+1.,su3vecAscendingComplex]
ForceFermionStaggeredOdd[gaugeMomOne,nonTrivial3x3mat,-1.,su3vecAscendingComplex]


(*force:Fermion Force Staggered with gaugeMomOne and IdentityMatrix[3]*)
(*Even, su3vecCold*)
ForceFermionStaggeredEven[gaugeMomOne,IdentityMatrix[3],1.,su3vecCold]
ForceFermionStaggeredEven[gaugeMomOne,IdentityMatrix[3],-1.,su3vecCold]
(*Odd, su3vecCold*)
ForceFermionStaggeredOdd[gaugeMomOne,IdentityMatrix[3],1.,su3vecCold]
ForceFermionStaggeredOdd[gaugeMomOne,IdentityMatrix[3],-1.,su3vecCold]
(*Even, su3vecAscendingComplex*)
ForceFermionStaggeredEven[gaugeMomOne,IdentityMatrix[3],+1.,su3vecAscendingComplex]
ForceFermionStaggeredEven[gaugeMomOne,IdentityMatrix[3],-1.,su3vecAscendingComplex]
(*Odd, su3vecAscendingComplex*)
ForceFermionStaggeredOdd[gaugeMomOne,IdentityMatrix[3],+1.,su3vecAscendingComplex]
ForceFermionStaggeredOdd[gaugeMomOne,IdentityMatrix[3],-1.,su3vecAscendingComplex]


(*force:Fermion Force Staggered with gaugeMomAscending and nonTrivial3x3mat*)
(*Even, su3vecCold*)
ForceFermionStaggeredEven[gaugeMomAscending,nonTrivial3x3mat,1.,su3vecCold]
ForceFermionStaggeredEven[gaugeMomAscending,nonTrivial3x3mat,-1.,su3vecCold]
(*Odd, su3vecCold*)
ForceFermionStaggeredOdd[gaugeMomAscending,nonTrivial3x3mat,1.,su3vecCold]
ForceFermionStaggeredOdd[gaugeMomAscending,nonTrivial3x3mat,-1.,su3vecCold]
(*Even, su3vecAscendingComplex*)
ForceFermionStaggeredEven[gaugeMomAscending,nonTrivial3x3mat,+1.,su3vecAscendingComplex]
ForceFermionStaggeredEven[gaugeMomAscending,nonTrivial3x3mat,-1.,su3vecAscendingComplex]
(*Odd, su3vecAscendingComplex*)
ForceFermionStaggeredOdd[gaugeMomAscending,nonTrivial3x3mat,+1.,su3vecAscendingComplex]
ForceFermionStaggeredOdd[gaugeMomAscending,nonTrivial3x3mat,-1.,su3vecAscendingComplex]


(*force:Fermion Force Staggered with gaugeMomOne and IdentityMatrix[3]*)
(*Even, su3vecCold*)
ForceFermionStaggeredEven[gaugeMomAscending,IdentityMatrix[3],1.,su3vecCold]
ForceFermionStaggeredEven[gaugeMomAscending,IdentityMatrix[3],-1.,su3vecCold]
(*Odd, su3vecCold*)
ForceFermionStaggeredOdd[gaugeMomAscending,IdentityMatrix[3],1.,su3vecCold]
ForceFermionStaggeredOdd[gaugeMomAscending,IdentityMatrix[3],-1.,su3vecCold]
(*Even, su3vecAscendingComplex*)
ForceFermionStaggeredEven[gaugeMomAscending,IdentityMatrix[3],+1.,su3vecAscendingComplex]
ForceFermionStaggeredEven[gaugeMomAscending,IdentityMatrix[3],-1.,su3vecAscendingComplex]
(*Odd, su3vecAscendingComplex*)
ForceFermionStaggeredOdd[gaugeMomAscending,IdentityMatrix[3],+1.,su3vecAscendingComplex]
ForceFermionStaggeredOdd[gaugeMomAscending,IdentityMatrix[3],-1.,su3vecAscendingComplex]



