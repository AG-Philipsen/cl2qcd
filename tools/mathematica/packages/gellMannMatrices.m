(* ::Package:: *)

BeginPackage["GellMannMatrices`"]

Lambda1::usage = 
	"Lambda1 gives the 1st Gell-Mann matrix."

Lambda2::usage = 
	"Lambda2 gives the 2nd Gell-Mann matrix."

Lambda3::usage = 
	"Lambda3 gives the 3rd Gell-Mann matrix."

Lambda4::usage = 
	"Lambda4 gives the 4th Gell-Mann matrix."

Lambda5::usage = 
	"Lambda5 gives the 5th Gell-Mann matrix."

Lambda6::usage = 
	"Lambda6 gives the 6th Gell-Mann matrix."

Lambda7::usage = 
	"Lambda7 gives the 7th Gell-Mann matrix."

Lambda8::usage = 
	"Lambda8 gives the 8th Gell-Mann matrix."

Begin["Private`"]

Lambda1:=
	Module[ {Lambda1={{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}},
	Lambda1
	]

Lambda2:=
	Module[ {Lambda2={{0, -I, 0}, {I, 0, 0}, {0, 0, 0}}},
	Lambda2
	]

Lambda3:=
	Module[ {Lambda3={{1, 0, 0}, {0, -1, 0}, {0, 0, 0}}},
	Lambda3
	]

Lambda4:=
	Module[ {Lambda4={{0, 0, 1}, {0, 0, 0}, {1, 0, 0}}},
	Lambda4
	]

Lambda5:=
	Module[ {Lambda5={{0, 0, -I}, {0, 0, 0}, {I, 0, 0}}},
	Lambda5
	]

Lambda6:=
	Module[ {Lambda6={{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}},
	Lambda6
	]

Lambda7:=
	Module[ {Lambda7={{0, 0, 0}, {0, 0, -I}, {0, I, 0}}},
	Lambda7
	]

Lambda8:=
	Module[ {Lambda8={{1, 0, 0}, {0, 1, 0}, {0, 0, -2}}/Sqrt[3]},
	Lambda8
	]

End[]

EndPackage[]
