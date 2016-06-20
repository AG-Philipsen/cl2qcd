(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]
Get["Packages.m",Path->{NotebookDirectory[]}]


(* ::Input:: *)
(*(*Fermionmatrices for a specific site,for constant gauge- and spinorfield*)Dslash[k_,s_,u_]:=k*((UnitMatrix[u]-Gamma1[u]).s+(UnitMatrix[u]-Gamma2[u]).s+(UnitMatrix[u]-Gamma3[u]).s+(UnitMatrix[u]-Gamma4[u]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma1[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma2[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma3[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma4[u\[ConjugateTranspose]]).s)*)
(**)
(*DslashBC[k_,s_,u_,ns_,nt_,ths_,tht_]:=k*((UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma1[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma2[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma3[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[tht*I*Pi/nt]]-Gamma4[u*Exp[tht*I*Pi/nt]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma1[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma2[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma3[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]+Gamma4[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]).s)*)
(**)
(*DslashBCAndChemPotIm[k_,s_,u_,ns_,nt_,ths_,tht_,cp_]:=k*((UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma1[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma2[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma3[u*Exp[ths*I*Pi/ns]]).s+Exp[+I*cp]*(UnitMatrix[u*Exp[tht*I*Pi/nt]]-Gamma4[u*Exp[tht*I*Pi/nt]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma1[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma2[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma3[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+Exp[+I*cp]*(UnitMatrix[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]+Gamma4[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]).s)*)
(**)
(*MWilson[k_,s_,u_]:=s-Dslash[k,s,u];*)
(**)
(*MTwistedMassPlus[k_,s_,u_,m_]:=(IdentityMatrix[12]+2. I*k*m*Gamma5[IdentityMatrix[3]]).s-Dslash[k,s,u];*)
(**)
(*MTwistedMassMinus[k_,s_,u_,m_]:=(IdentityMatrix[12]-2.*I*k*m*Gamma5[IdentityMatrix[3]]).s-Dslash[k,s,u];*)
(**)
(*MTwistedMassSitediagonal[s_,k_,m_]:=(IdentityMatrix[12]+2.*I*m*k*Gamma5[IdentityMatrix[3]]).s;*)
(**)
(*MTwistedMassSitediagonalMinus[s_,k_,m_]:=(IdentityMatrix[12]-2.*I*m*k*Gamma5[IdentityMatrix[3]]).s;*)
(**)
(*MTwistedMassInverseSitediagonal[s_,k_,m_]:=(1./(1.+4*k*k*m*m))*(IdentityMatrix[12]-2.*I*m*k*Gamma5[IdentityMatrix[3]]).s;*)
(**)
(*MTwistedMassInverseSitediagonalMinus[s_,k_,m_]:=(1./(1.+4*k*k*m*m))*(IdentityMatrix[12]+2.*I*m*k*Gamma5[IdentityMatrix[3]]).s;*)
(**)
(*SaxpyAndGamma5EvenOdd[s1_,s2_,a_]:=Gamma5[IdentityMatrix[3]].(s2-a*s1);*)
(**)
(*MTwistedMassSitediagonalAndGamma5EvenOdd[s_,k_,m_]:=Gamma5[IdentityMatrix[3]].MTwistedMassSitediagonal[s,k,m]*)
(**)
(*MTwistedMassSitediagonalMinusAndGamma5EvenOdd[s_,k_,m_]:=Gamma5[IdentityMatrix[3]].MTwistedMassSitediagonalMinus[s,k,m]*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalMinus[s_,u_,k_,m_]:=MTwistedMassInverseSitediagonalMinus[(-1.)*Dslash[k,s,u],k,m]*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonal[s_,u_,k_,m_]:=MTwistedMassInverseSitediagonal[(-1.)*Dslash[k,s,u],k,m]*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalBC[s_,u_,k_,m_,ns_,nt_,ths_,tht_]:=MTwistedMassInverseSitediagonal[(-1.)*DslashBC[k,s,u,ns,nt,ths,tht],k,m]*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[s_,u_,k_,m_,ns_,nt_,ths_,tht_,cp_]:=MTwistedMassInverseSitediagonal[(-1.)*DslashBCAndChemPotIm[k,s,u,ns,nt,ths,tht,cp],k,m]*)
(**)
(*DslashBCAndChemPotImAndRe[k_,s_,u_,ns_,nt_,ths_,tht_,cpRe_,cpIm_]:=k*((UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma1[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma2[u*Exp[ths*I*Pi/ns]]).s+(UnitMatrix[u*Exp[ths*I*Pi/ns]]-Gamma3[u*Exp[ths*I*Pi/ns]]).s+Exp[+I*cpIm]*Exp[+cpRe]*(UnitMatrix[u*Exp[tht*I*Pi/nt]]-Gamma4[u*Exp[tht*I*Pi/nt]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma1[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma2[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+(UnitMatrix[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]+Gamma3[u\[ConjugateTranspose]*Exp[-ths*I*Pi/ns]]).s+Exp[-I*cpIm]*Exp[-cpRe]*(UnitMatrix[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]+Gamma4[u\[ConjugateTranspose]*Exp[-tht*I*Pi/nt]]).s)*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[s_,u_,k_,m_,ns_,nt_,ths_,tht_,cpIm_,cpRe_]:=MTwistedMassInverseSitediagonal[(-1.)*DslashBCAndChemPotImAndRe[k,s,u,ns,nt,ths,tht,cpIm,cpRe],k,m]*)
(**)
(*DslashChemPotIm[k_,s_,u_,cpIm_]:=k*((UnitMatrix[u]-Gamma1[u]).s+(UnitMatrix[u]-Gamma2[u]).s+(UnitMatrix[u]-Gamma3[u]).s+Exp[+I*cpIm]*(UnitMatrix[u]-Gamma4[u]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma1[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma2[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma3[u\[ConjugateTranspose]]).s+Exp[-I*cpIm]*(UnitMatrix[u\[ConjugateTranspose]]+Gamma4[u\[ConjugateTranspose]]).s)*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalChemPotIm[s_,u_,k_,m_,ns_,nt_,cpIm_]:=MTwistedMassInverseSitediagonal[(-1.)*DslashChemPotIm[k,s,u,cpIm],k,m]*)
(**)
(*DslashChemPotRe[k_,s_,u_,cpRe_]:=k*((UnitMatrix[u]-Gamma1[u]).s+(UnitMatrix[u]-Gamma2[u]).s+(UnitMatrix[u]-Gamma3[u]).s+Exp[cpRe]*(UnitMatrix[u]-Gamma4[u]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma1[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma2[u\[ConjugateTranspose]]).s+(UnitMatrix[u\[ConjugateTranspose]]+Gamma3[u\[ConjugateTranspose]]).s+Exp[cpRe]*(UnitMatrix[u\[ConjugateTranspose]]+Gamma4[u\[ConjugateTranspose]]).s)*)
(**)
(*DslashEvenOddAndMTmInverseSitediagonalChemPotRe[s_,u_,k_,m_,ns_,nt_,cpRe_]:=MTwistedMassInverseSitediagonal[(-1.)*DslashChemPotRe[k,s,u,cpRe],k,m]*)
(**)


(*Tests for MWilson checking the squarenorm of the resulting spinor*)
squareNorm[MWilson[0., spinorCold, coldSU3mat]]
squareNorm[MWilson[0.,spinorAscendingComplex,coldSU3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorCold,coldSU3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorAscendingComplex,coldSU3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorCold,nonTrivialSU3mat]]
squareNorm[MWilson[nonTrivialRealPar,spinorAscendingComplex,nonTrivialSU3mat]]


(*Tests for MWilson checking the sum of the resulting spinor*)
countSf[MWilson[0.,spinorCold,coldSU3mat]]
countSf[MWilson[0.,spinorAscendingComplex,coldSU3mat]]
countSf[MWilson[nonTrivialRealPar,spinorCold,coldSU3mat]]
countSf[MWilson[nonTrivialRealPar,spinorAscendingComplex,coldSU3mat]]
countSf[MWilson[nonTrivialRealPar,spinorCold,nonTrivialSU3mat]]
countSf[MWilson[nonTrivialRealPar,spinorAscendingComplex,nonTrivialSU3mat]]


(*Tests for MTwistedMassMinus checking the squarenorm of the resulting spinor*)
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,coldSU3mat,0.]]
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,coldSU3mat,nonTrivialRealPar]]
squareNorm[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar]]


(*Tests for MTwistedMassMinus checking the sum of the resulting spinor*)
countSf[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,coldSU3mat,nonTrivialRealPar]]
countSf[MTwistedMassMinus[nonTrivialRealPar,spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar]]


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
countSf[DslashEvenOddAndMTmInverseSitediagonal[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,0.]]
countSf[DslashEvenOddAndMTmInverseSitediagonal[spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar,0.]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBC checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBC[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBC[spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPot[spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and boundary condition phases (ths,tht)*)
FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar,nonTrivialRealPar,ns,nt,ths,tht,nonTrivialRealPar,nonTrivialRealPar]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals}]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt), the boundary condition phases (ths,tht), and the chemical potential (cpRe,cpIm)*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]
(*Imaginary part of the result:*)
ComplexExpand[Im[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalChemPotIm checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and the imaginary chemical potential cpIm*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalChemPotIm[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,0.,ns,nt,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,cpIm\[Element]Reals}]]]]
(*Tests for DslashEvenOddAndMTmInverseSitediagonalChemPotRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt) and the real chemical potential cpRe*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalChemPotRe[spinorAscendingComplex,coldSU3mat,nonTrivialRealPar,0.,ns,nt,cpRe]],{ns\[Element]Reals,nt\[Element]Reals,cpRe\[Element]Reals}]]]]


(*Tests for DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe checking the sum of the resulting spinor and expressing it as a function of the lattice extents (ns, nt), the boundary condition phases (ths,tht) and the chemical potential (cpRe,cpIm)*)
(*Real part of the result:*)
ComplexExpand[Re[ExpToTrig[FullSimplify[countSf[DslashEvenOddAndMTmInverseSitediagonalBCAndChemPotImAndRe[spinorAscendingComplex,nonTrivialSU3mat,nonTrivialRealPar,0.,ns,nt,ths,tht,cpRe,cpIm]],{ns\[Element]Reals,nt\[Element]Reals,ths\[Element]Reals,tht\[Element]Reals,cpRe\[Element]Reals,cpIm\[Element]Reals}]]]]
