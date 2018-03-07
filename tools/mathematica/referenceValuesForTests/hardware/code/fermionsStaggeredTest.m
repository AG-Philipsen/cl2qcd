(* ::Package:: *)

(*
 * Copyright 2012,2013 Lars Zeidlewicz,Christopher Pinke,
 * Matthias Bach,Christian Sch\[ADoubleDot]fer,Stefano Lottini,Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
*)


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
