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


stripMinus[x_]:=x/.{-a_->a};
replaceCD[expr_]:=expr/.{PD[idx_][val_]:>CD[stripMinus[idx],replaceCD[val]],CD[idx_][val_]:>CD[stripMinus[idx],replaceCD[val]]};
printAndReturn[expr_]:=(Print[ToString[expr]];expr);
LI=xAct`xTensor`LI;
XActToMatcheteInternal[expr_,translationRules_]:=replaceCD[(expr /. {hh_[LI[order_],idx1_,idx2_]
															:>If[order==1,(hh/.translationRules)[stripMinus[idx1],stripMinus[idx2]],0]}/.{
															hh_[idx1_,idx2_]/;MatchQ[hh, Alternatives @@ (First /@ translationRules)]:>
															(hh/.translationRules)[stripMinus[idx1],stripMinus[idx2]]})];

XActToMatchete[expr_Plus,translationRules_] := Plus @@ ((XActToMatchete[#,translationRules] &) /@ List @@ expr);
XActToMatchete[expr_,translationRules_] := Module[{e = expr},
  e = translateHeads[e, xActToMatcheteHeads];
  XActToMatcheteInternal[e,translationRules]
];

End[];
EndPackage[];
