(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]]], "packages/threeBythreeMatrices.m"}], Path -> {NotebookDirectory[]}]


(*Staple_sum*)
stSum[u_] := 
  3*(u.ConjugateTranspose[u].ConjugateTranspose[u] + 
        ConjugateTranspose[u].ConjugateTranspose[u].u)
(*Calc_staple tests*)
(Re[Total[stSum[cold3x3mat], {1, 2}]] + Im[Total[stSum[cold3x3mat], {1, 2}]]) * 4
(Re[Total[stSum[nonTrivial3x3mat], {1, 2}]] + Im[Total[stSum[nonTrivial3x3mat], {1, 2}]]) * 4



