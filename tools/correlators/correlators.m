(* ::Package:: *)

(* ::Input:: *)
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
