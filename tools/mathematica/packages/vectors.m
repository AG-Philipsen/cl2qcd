(* ::Package:: *)

(*
 * Copyright (c) 2016 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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


BeginPackage["Vectors`"]

squareNorm::usage =
	"squareNorm[v] takes a generic vector as argument and compute its norm squared."

Begin["Private`"]

squareNorm[v_]:=
	Module[{squareNorm=Norm[v]^2},
	squareNorm
	]

End[]

EndPackage[]
