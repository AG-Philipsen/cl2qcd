(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/staggeredDiracOperator.m"}], Path -> {NotebookDirectory[]}]


Chop[FullSimplify[
  MStaggeredBC[cold3x3mat, su3vecAscendingComplex, 
   nonTrivialRealPar, tht, ths, nt, ns], {tht \[Element] Reals, 
   nt \[Element] Reals, ths \[Element] Reals, ns \[Element] Reals}], 
 10^-10]


Chop[ExpToTrig[
  FullSimplify[
   MStaggeredBC[nonTrivial3x3mat, su3vecCold, nonTrivialRealPar,
     tht, ths, nt, ns], {tht \[Element] Reals, nt \[Element] Reals, 
    ths \[Element] Reals, ns \[Element] Reals}]], 10^-10]


Chop[FullSimplify[
  MStaggeredBC[nonTrivial3x3mat, su3vecAscendingComplex, nonTrivialRealPar,
	tht, ths, nt, ns], {tht \[Element] Reals, nt \[Element] Reals,
	ths \[Element] Reals, ns \[Element] Reals}], 10^-10]


Chop[FullSimplify[
  DKStaggered[cold3x3mat, su3vecAscendingComplex, 0., cp, tht, ths,
    nt, ns], {cp \[Element] Reals, tht \[Element] Reals, 
   nt \[Element] Reals, ths \[Element] Reals, ns \[Element] Reals}], 
 10^-10]


Chop[ExpToTrig[
  FullSimplify[
   DKStaggered[nonTrivial3x3mat, su3vecCold, 0., cp, tht, ths, nt, 
    ns], {cp \[Element] Reals, tht \[Element] Reals, 
    nt \[Element] Reals, ths \[Element] Reals, ns \[Element] Reals}]],
  10^-10]


Chop[FullSimplify[
  DKStaggered[nonTrivial3x3mat, su3vecAscendingComplex, 0., cp, tht, ths,
    nt, ns], {cp \[Element] Reals, tht \[Element] Reals, 
   nt \[Element] Reals, ths \[Element] Reals, ns \[Element] Reals}], 
 10^-10]
