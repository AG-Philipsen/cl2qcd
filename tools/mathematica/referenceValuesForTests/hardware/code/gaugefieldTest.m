(* ::Package:: *)

SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]], "packages/threeBythreeMatrices.m"}], Path -> {NotebookDirectory[]}]


(*Plaquette*)
Plaquette[u_] := Tr[u.u.ConjugateTranspose[u].ConjugateTranspose[u]]

(*Plaquette tests*)
Plaquette[cold3x3mat]
Plaquette[nonTrivial3x3mat]


(*Rectangles*)
Rect[u_] := 
 Tr[u.u.u.ConjugateTranspose[u].ConjugateTranspose[
      u].ConjugateTranspose[u]]/3*6*2

(*Rectangles test*)
Rect[cold3x3mat]
Rect[nonTrivial3x3mat]


(*Polyakov for nt4*)
Polyakov[u_] := Tr[u.u.u.u]/3

(*Polyakov test*)
Polyakov[cold3x3mat]
Polyakov[nonTrivial3x3mat]



