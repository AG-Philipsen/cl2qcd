(* ::Package:: *)

Get["gellMannMatrices.m", Path -> {NotebookDirectory[]}]
BeginPackage["EightComponentsVectors`"]

gaugeMomZero::usage = 
	"gaugeMomZero gives a zero algebra vector (uniformly filled by {0.})."

gaugeMomOne::usage = 
	"gaugeMomOne gives a cold algebra vector (uniformly filled by {1.})."

gaugeMomAscending::usage =
	"gaugeMomAscending gives an algebra vector filled with ascending numbers."

algebraElement::usage =
	"algebraElement takes an SU(3) matrix as argument and computes the algebra element
	 Tr[i * T_k * (in - in^\\dag)], i.e. calculates the trace of I times a generator of the
	 algebra su(3) times a 3x3-matrix (see tr_lambda_u in operations_gaugemomentum.cl)."

aeFromSu3::usage = 
	"aeFromSu3 takes a 3-by-3 matrix and builds an algebra element out of it
	(see build_ae_from_su3 in operations_matrix_su3.cl)."

Begin["Private`"]
Needs["GellMannMatrices`"]

gaugeMomZero:=
	Module[ {gaugeMomZero={0., 0., 0., 0., 0., 0., 0., 0.}},
	gaugeMomZero
	]

gaugeMomOne:=
	Module[ {gaugeMomOne={1., 1., 1., 1., 1., 1., 1., 1.}},
	gaugeMomOne
	]

gaugeMomAscending:=
	Module[ {gaugeMomAsc={1., 2., 3., 4., 5., 6., 7., 8.}},
	gaugeMomAsc
	]

algebraElement[u_]:=
	Module[{algebraElement={Tr[I*Lambda1.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda2.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda3.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda4.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda5.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda6.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda7.(u - ConjugateTranspose[u])]/2,
							Tr[I*Lambda8.(u - ConjugateTranspose[u])]/2}},
	algebraElement
	]

aeFromSu3[u_]:=
	Module[{aeFromSu3={2*Re[u[[1, 2]]], -2*Im[u[[1, 2]]], Re[u[[1, 1]]] - Re[u[[2, 2]]],
						2*Re[u[[1, 3]]], -2 Im[u[[1, 3]]], 2*Re[u[[2, 3]]], -2*Im[u[[2, 3]]],
						(Re[u[[1, 1]]] + Re[u[[2, 2]]])*3/Sqrt[3]}},
	aeFromSu3
	]

End[]

EndPackage[]

