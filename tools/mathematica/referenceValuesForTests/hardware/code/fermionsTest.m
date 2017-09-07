(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/wilsonDiracOperator.m"}], Path -> {NotebookDirectory[]}]


(*Tests for MWilson checking the squarenorm of the resulting spinor*)
squareNorm[MWilson[0., spinorCold, cold3x3mat]]
squareNorm[MWilson[0.,spinorAscendingComplex,cold3x3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorCold,cold3x3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorAscendingComplex,cold3x3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorCold,nonTrivial3x3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat]]


(*Tests for MWilson checking the sum of the resulting spinor*)
countSf[MWilson[0.,spinorCold,cold3x3mat]]
countSf[MWilson[0.,spinorAscendingComplex,cold3x3mat]]
countSf[MWilson[nonTrivialRealPar,spinorCold,cold3x3mat]]
countSf[MWilson[nonTrivialRealPar,spinorAscendingComplex,cold3x3mat]]
countSf[MWilson[nonTrivialRealPar,spinorCold,nonTrivial3x3mat]]
countSf[MWilson[nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat]]


(*Tests for MTwistedMassMinus checking the squarenorm of the resulting spinor*)
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,cold3x3mat,0.]]
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,cold3x3mat,nonTrivialRealPar]]
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar]]


(*Tests for MTwistedMassMinus checking the sum of the resulting spinor*)
countSf[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,cold3x3mat,nonTrivialRealPar]]
countSf[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar]]


(*Tests for MTwistedMassSitediagonal checking the squarenorm of the resulting spinor*)
squareNorm[MTwistedMassSitediagonal[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]
squareNorm[MTwistedMassInverseSitediagonal[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]


(*Tests for MTwistedMassSitediagonal checking the sum of the resulting spinor*)
countSf[MTwistedMassSitediagonal[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]
countSf[MTwistedMassInverseSitediagonal[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]


(*Tests for MTwistedMassSitediagonalMinus checking the squarenorm of the resulting spinor*)
squareNorm[MTwistedMassSitediagonalMinus[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]
squareNorm[MTwistedMassInverseSitediagonalMinus[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]


(*Tests for MTwistedMassSitediagonalMinus checking the sum of the resulting spinor*)
countSf[MTwistedMassSitediagonalMinus[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]
countSf[MTwistedMassInverseSitediagonalMinus[spinorAscendingComplex,nonTrivialRealPar,nonTrivialRealPar]]


(*Tests for SaxpyAndGamma5EvenOdd checking the sum of the resulting spinor*)
countSf[SaxpyAndGamma5EvenOdd[spinorAscendingComplex,spinorAscendingComplex,0.+I*0.]]
countSf[SaxpyAndGamma5EvenOdd[spinorAscendingComplex,spinorAscendingComplex,(-1.)*nonTrivialRealPar+I*2.*nonTrivialRealPar]]
countSf[Gamma5[IdentityMatrix[3]].spinorAscendingComplex]


(*Tests for MTwistedMassSitediagonalAndGamma5EvenOdd checking the sum of the resulting spinor*)
countSf[MTwistedMassSitediagonalAndGamma5EvenOdd[spinorAscendingComplex,0.,0.]]
countSf[MTwistedMassSitediagonalAndGamma5EvenOdd[spinorAscendingComplex,-4.*nonTrivialRealPar,11*nonTrivialRealPar]]
countSf[MTwistedMassSitediagonalMinusAndGamma5EvenOdd[spinorAscendingComplex,0.,0.]]
countSf[MTwistedMassSitediagonalMinusAndGamma5EvenOdd[spinorAscendingComplex,2.*nonTrivialRealPar,-5*nonTrivialRealPar]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonal checking the sum of the resulting spinor*)
countSf[DslashEvenOddAndMTmInverseSitediagonal[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,0.]]
countSf[DslashEvenOddAndMTmInverseSitediagonal[spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar,0.]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBC checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBC[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBC[spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt), the boundary condition phases (ths,tht), and the chemical potential (cpRe,cpIm)*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]
(*Imaginary part of the result:*)
ComplexExpand[Im[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalChemPotIm checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and the imaginary chemical potential cpIm*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalChemPotIm[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,0.,ns,nt,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,cpIm\[Element]Reals}]]]]
(*Tests for DslashEvenOddAndMTmInverseSitediagonalChemPotRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and the real chemical potential cpRe*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalChemPotRe[spinorAscendingComplex,cold3x3mat,nonTrivialRealPar,0.,ns,nt,cpRe]],{ns\[Element]Reals,nt\[Element]Reals,cpRe\[Element]Reals}]]]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt), the boundary condition phases (ths,tht) and the chemical potential (cpRe,cpIm)*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,nonTrivial3x3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]
