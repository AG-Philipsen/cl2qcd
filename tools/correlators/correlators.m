(* ::Package:: *)

(* ::Input:: *)
(*(**)
(* * Copyright 2012,2013 Lars Zeidlewicz,Christopher Pinke,*)
(* * Matthias Bach,Christian Sch\[ADoubleDot]fer,Stefano Lottini,Alessandro Sciarra*)
(* **)
(* * This file is part of CL2QCD.*)
(* **)
(* * CL2QCD is free software:you can redistribute it and/or modify*)
(* * it under the terms of the GNU General Public License as published by*)
(* * the Free Software Foundation,either version 3 of the License,or*)
(* * (at your option) any later version.*)
(* **)
(* * CL2QCD is distributed in the hope that it will be useful,*)
(* * but WITHOUT ANY WARRANTY;without even the implied warranty of*)
(* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the*)
(* * GNU General Public License for more details.*)
(* **)
(* * You should have received a copy of the GNU General Public License*)
(* * along with CL2QCD.If not,see<http://www.gnu.org/licenses/>.*)
(**)*)
(**)
(*gamma5:=({*)
(* {1, 0, 0, 0},*)
(* {0, 1, 0, 0},*)
(* {0, 0, -1, 0},*)
(* {0, 0, 0, -1}*)
(*})*)


(* ::Input:: *)
(*gamma4:=({*)
(* {0, 0, -1, 0},*)
(* {0, 0, 0, -1},*)
(* {-1, 0, 0, 0},*)
(* {0, -1, 0, 0}*)
(*})*)


(* ::Input:: *)
(*gamma3:=({*)
(* {0, 0, -I, 0},*)
(* {0, 0, 0, I},*)
(* {I, 0, 0, 0},*)
(* {0, -I, 0, 0}*)
(*})*)


(* ::Input:: *)
(*gamma2:=({*)
(* {0, 0, 0, -1},*)
(* {0, 0, 1, 0},*)
(* {0, 1, 0, 0},*)
(* {-1, 0, 0, 0}*)
(*})*)


(* ::Input:: *)
(*gamma1:=({*)
(* {0, 0, 0, -I},*)
(* {0, 0, -I, 0},*)
(* {0, I, 0, 0},*)
(* {I, 0, 0, 0}*)
(*})*)


(* ::Input:: *)
(*unit:=({*)
(* {1, 0, 0, 0},*)
(* {0, 1, 0, 0},*)
(* {0, 0, 1, 0},*)
(* {0, 0, 0, 1}*)
(*})*)


(* ::Input:: *)
(*anticomm[A_,B_]:=A.B+B.A*)


(* ::Input:: *)
(*correlator[G_,v_]:=-Tr[(Conjugate[v].G).(Transpose[v].G)]*)


(* ::Input:: *)
(*gamma41:=gamma4.gamma1*)


(* ::Input:: *)
(*gamma42:=gamma4.gamma2*)


(* ::Input:: *)
(*gamma43:=gamma4.gamma3*)


(* ::Input:: *)
(*gamma541:=gamma5.gamma41*)


(* ::Input:: *)
(*gamma542:=gamma5.gamma42*)


(* ::Input:: *)
(*gamma543:=gamma5.gamma43*)


(* ::Input:: *)
(*phi:=({*)
(* {p00, p01, p02, p03},*)
(* {p10, p11, p12, p13},*)
(* {p20, p21, p22, p23},*)
(* {p30, p31, p32, p33}*)
(*})*)


(* ::Input:: *)
(*correlator[unit,phi]*)


(* ::Input:: *)
(*correlator[gamma5,phi]*)


(* ::Input:: *)
(*correlator[gamma43,phi]*)
