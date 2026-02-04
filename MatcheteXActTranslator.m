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
PD=xAct`xTensor`PD;
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

sanitizeIndices[idxs_]:=If[MatchQ[idxs, {idx_,idx_}],{idxs[[1]],-idxs[[2]]},idxs]

(* --- Handle individual Matchete constructs --- *)
translateField[Field[name_, type_, {}, derivsRaw_List],translationRules_] := Module[
  {idxsRaw, derivs, idxs, translatedName, base},
  (* Extract only the symbol part from Index[...] *)
  derivs = First /@ derivsRaw;
  (* For tensor fields, extract head arguments likewise *)
  idxsRaw = If[Head[type] === Graviton, List @@ type, {}];
  idxs = First /@ idxsRaw;
  (* Build the base tensor/scalar *)
  translatedName=name/.translationRules;
  base = If[idxs === {}, translatedName[], translatedName @@ sanitizeIndices[idxs]];
  (* Wrap partial derivatives *)
  Fold[Function[{expr, idx}, PD[idx][expr]], base, derivs]
];

translateCoupling[Coupling[name_, {}, _]] := name;

lowerAllIndices[expr_] := expr /. {PD[i_][sub_]:>PD[-i][lowerAllIndices[sub]],h_Symbol[idx1_,idx2_]:>h[-idx1,-idx2]}

(* --- Avoid index clashes by alternating signs --- *)
sanitizeMonomial[term_] := Module[{counts, fIdx, excludeHeads,indexSubstitution},
  (* Initialize counter for indices and list of non-index heads *)
  counts = <||>;
  excludeHeads = {Plus, Times, PD, Rational, Power};

  (* fIdx: on first use returns original index, afterwards returns negative to denote lowered index *)
  fIdx[idx_Symbol] := (
    counts[idx] = Lookup[counts, idx, 0] + 1;
    (*Print["index ",idx," counts",counts[idx]];*)
    If[counts[idx] == 1, idx, -idx]
  );

  (* Perform substitutions: flip in PD wrappers and tensor heads *)
  indexSubstitution[term1_] := (
  term1 /. {
    (* For partial derivatives PD[index][expr] *)
    PD[i_][sub_] :> (PD[fIdx[i]])[indexSubstitution[sub]],

    (* For any tensor head h[arg1,arg2,...] where all args are symbols and head not excluded *)
    h_Symbol[inds__] /; (! MemberQ[excludeHeads, h] && And @@ (MatchQ[#, _Symbol] & /@ {inds})) :>
      h @@ (fIdx /@ {inds})
  }
  );
  indexSubstitution[term] /.{h_Symbol[idx1_,idx2_]^2:>h[idx1,idx2]h[-idx1,-idx2],PD[i_][sub_]^2:>PD[i][sub]lowerAllIndices[PD[i][sub]]}
];

(* --- Main Matchete -> xAct translator --- *)
translateTermMatchete[term_,translationRules_] := Module[{e = term},
  e = e //. {f_Field :> translateField[f,translationRules], c_Coupling :> translateCoupling[c]};
  translateHeads[e, matcheteToXActHeads]
];

MatcheteToXActInternal[expr_Plus,translationRules_] := Plus @@ (sanitizeMonomial[translateTermMatchete[#,translationRules]] & /@ List @@ expr);
MatcheteToXActInternal[expr_,translationRules_]      := sanitizeMonomial[translateTermMatchete[expr,translationRules]];
MatcheteToXAct[expr_,translationRules_] := MatcheteToXActInternal[Distribute[expr, Plus], translationRules];

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
      Field[base, Matchete`Scalar, {}, derivs],
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
