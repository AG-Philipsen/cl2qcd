(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]]], "packages/threeBythreeMatrices.m"}], Path -> {NotebookDirectory[]}]


(*local_Q*)
localQ[u_] := Tr[u.u.ConjugateTranspose[u].ConjugateTranspose[u]]*6*4

(*local_Q tests(see miscellaneousTest)*)
localQ[cold3x3mat]
localQ[nonTrivial3x3mat]



