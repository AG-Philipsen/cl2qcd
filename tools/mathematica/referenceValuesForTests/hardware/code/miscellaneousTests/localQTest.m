(* ::Package:: *)

(*
 * Copyright 2012,2013 Lars Zeidlewicz,Christopher Pinke,
 * Matthias Bach,Christian Sch\[ADoubleDot]fer,Stefano Lottini,Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software:you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation,either version 3 of the License,or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY;without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.If not,see<http://www.gnu.org/licenses/>.
*)


SetOptions[SelectedNotebook[], PrintPrecision -> 16]

(*Get the needed packages for these tests*)
Get[FileNameJoin[{ParentDirectory[ParentDirectory[ParentDirectory[ParentDirectory[NotebookDirectory[]]]]], "packages/threeBythreeMatrices.m"}], Path -> {NotebookDirectory[]}]


(*local_Q*)
localQ[u_] := Tr[u.u.ConjugateTranspose[u].ConjugateTranspose[u]]*6*4

(*local_Q tests(see miscellaneousTest)*)
localQ[cold3x3mat]
localQ[nonTrivial3x3mat]
