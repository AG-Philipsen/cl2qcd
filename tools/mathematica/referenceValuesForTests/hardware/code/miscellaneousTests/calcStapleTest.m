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


(*Staple_sum*)
stSum[u_] :=
  3*(u.ConjugateTranspose[u].ConjugateTranspose[u] +
        ConjugateTranspose[u].ConjugateTranspose[u].u)
(*Calc_staple tests*)
(Re[Total[stSum[cold3x3mat], {1, 2}]] + Im[Total[stSum[cold3x3mat], {1, 2}]]) * 4
(Re[Total[stSum[nonTrivial3x3mat], {1, 2}]] + Im[Total[stSum[nonTrivial3x3mat], {1, 2}]]) * 4
