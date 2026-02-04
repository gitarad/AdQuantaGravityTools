(* ::Package:: *)

(* ::Package:: *)

(* :Title: MatcheteXActTranslator *)
(* :Context: MatcheteXActTranslator` *)
(* :Author: ChatGPT *)
(* :Summary: Utilities to translate expressions between Matchete and xAct xTensor syntax,
   ensuring no index clashes by alternating upper/lower index positions.*)
(* :Package Version: 1.5 *)
(* :Mathematica Version: 12.0+ *)

BeginPackage["MatcheteXActTranslator`", {"Matchete`"}];

MatcheteToXAct::usage = 
  "MatcheteToXAct[expr] converts a Matchete expression to xAct xTensor syntax,\
   automatically alternating index signs to avoid duplicate upper/lower clashes.";
XActToMatchete::usage = 
  "XActToMatchete[expr] converts an xAct xTensor expression back to Matchete syntax.";

Begin["`Private`"];

(* --- HEAD MAPPINGS --- *)
matcheteToXActHeads = <|
  "LeviCivitaSymbol" -> "EpsilonSymbol",
  "Metric"           -> "Metric",
  "KroneckerDelta"   -> "Delta",
  "CovariantD"       -> "CD",
  "PartialD"         -> "PD",
  "ChristoffelSymbol"-> "Christoffel",
  "RicciTensor"      -> "Ricci",
  "RiemannTensor"    -> "Riemann"
|>;
xActToMatcheteHeads = Association[Reverse /@ Normal[matcheteToXActHeads]];

translateHeads[expr_, mapping_] := expr /. {
    h_Symbol[args__] /; KeyExistsQ[mapping, SymbolName[h]] :>
      ToExpression[mapping[SymbolName[h]]] @@ (translateHeads[# , mapping] & /@ {args})
};

(* --- Handle individual Matchete constructs --- *)
translateField[Field[name_, type_, {}, derivsRaw_List]] := Module[
  {idxsRaw, derivs, idxs, base},
  (* Extract only the symbol part from Index[...] *)
  derivs = First /@ derivsRaw;
  (* For tensor fields, extract head arguments likewise *)
  idxsRaw = If[Head[type] === Graviton, List @@ type, {}];
  idxs = First /@ idxsRaw;
  (* Build the base tensor/scalar *)
  base = If[idxs === {}, name[], name @@ idxs];
  (* Wrap partial derivatives *)
  Fold[Function[{expr, idx}, PD[idx][expr]], base, derivs]
];

translateCoupling[Coupling[name_, {}, _]] := name;

(* --- Avoid index clashes by alternating signs --- *)
sanitizeMonomial[term_] := Module[{counts, fIdx, excludeHeads},
  (* Initialize counter for indices and list of non-index heads *)
  counts = <||>;
  excludeHeads = {Plus, Times, PD, Rational, Power};

  (* fIdx: on first use returns original index, afterwards returns negative to denote lowered index *)
  fIdx[idx_Symbol] := (
    counts[idx] = Lookup[counts, idx, 0] + 1;
    If[counts[idx] == 1, idx, -idx]
  );

  (* Perform substitutions: flip in PD wrappers and tensor heads *)
  term /. {
    (* For partial derivatives PD[index][expr] *)
    PD[i_][sub_] :> (PD[fIdx[i]])[sub],

    (* For any tensor head h[arg1,arg2,...] where all args are symbols and head not excluded *)
    h_Symbol[inds__] /; (! MemberQ[excludeHeads, h] && And @@ (MatchQ[#, _Symbol] & /@ {inds})) :>
      h @@ (fIdx /@ {inds})
  }
];

(* --- Main Matchete -> xAct translator --- *)
translateTermMatchete[term_] := Module[{e = term},
  e = e //. {f_Field :> translateField[f], c_Coupling :> translateCoupling[c]};
  translateHeads[e, matcheteToXActHeads]
];

MatcheteToXAct[expr_Plus] := Plus @@ (sanitizeMonomial[translateTermMatchete[#]] & /@ List @@ expr);
MatcheteToXAct[expr_]      := sanitizeMonomial[translateTermMatchete[expr]];

(* --- xAct -> Matchete translator remains unchanged --- *)
clearPD[expr_] := Module[{tmp = expr, d = {}},
  While[MatchQ[tmp, PD[_][_] ],
    d = Append[d, First@Replace[tmp, PD[i_][_] :> i]];
    tmp = Replace[tmp, PD[_][inner_] :> inner];
  ];
  {tmp, d}
];

translateExprXAct[expr_] := Module[{base, derivs},
  {base, derivs} = clearPD[expr];
  Which[
    Head[base] === Symbol,
      Field[base, Scalar, {}, derivs],
    True,
      Field[Head[base],
        If[Length[List @@ base] == 2,
           Graviton @@ (List @@ base),
           Head[base] @@ (List @@ base)
        ], {}, derivs]
  ]
];

XActToMatchete[expr_Plus] := Plus @@ (XActToMatchete /@ List @@ expr);
XActToMatchete[expr_] := Module[{e = expr},
  e = translateHeads[e, xActToMatcheteHeads];
  translateExprXAct[e]
];

End[];
EndPackage[];
