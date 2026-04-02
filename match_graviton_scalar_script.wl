(* ::Package:: *)

Print["Adding to path"]
PrependTo[$Path, "/home/gil.arad"]
Print["Loading Matchete"]
Needs["Matchete`"];
Print["Loading xAct"]
Needs["xAct`xTensor`"];
Needs["xAct`xPert`"];
Needs["xAct`xTras`"];
Needs["MatcheteXActTranslator`"];
Print["Defining field"]
DefineField[\[Phi], Matchete`Scalar, SelfConjugate -> True, 
  Mass -> {Heavy, M}];
DefineField[h, Graviton, SelfConjugate -> True, Mass -> {Light, 0}];
DefineCoupling[\[Kappa], SelfConjugate -> True];
Print[GetFields[]];
DefManifold[Mpert,4,{\[Alpha],\[Beta],\[Rho],\[Tau]}];
DefMetric[{1,3,0},g[-\[Alpha],-\[Beta]],CD,{";","\[Del]"}];
DefMetricPerturbation[g,H,\[CurlyEpsilon]];
DefTensor[\[CapitalPhi][],Mpert];
DefTensorPerturbation[\[CapitalPhi]pert[LI[0]],\[CapitalPhi][],\[CurlyEpsilon]];
DefManifold[M4,4,{a,b,c,d}];
DefMetric[{1,3,0},\[Eta][-a,-b],PD,{",","\[PartialD]"},FlatMetric->True];
DefTensor[r[-a,-b],M4,Symmetric[{-a,-b}]];
DefTensor[\[CurlyPhi][],M4];
DefConstantSymbol[{hbar,hbar2,\[Kappa],\[Kappa]2,\[Epsilon],M,\[Mu]bar2}];
perturbationOrder=7
Print["Expanding perturbation, order "<>ToString[perturbationOrder]]
perturbationSum[expr_,order_,orderParam_]:=Sum[PerturbFlat[expr,i]orderParam^i/Factorial[i],{i,1,order}]
perturbation = (ExpandPerturbation[perturbationSum[Sqrt[-Detg[]](g[\[Alpha],\[Beta]]CD[-\[Alpha]][\[CapitalPhi][]]CD[-\[Beta]][\[CapitalPhi][]]-M^2\[CapitalPhi][]^2)/2,perturbationOrder,\[Kappa]]]//ToFlat//CollectTensors)/.{Sqrt[-Detg[]]->1,\[CapitalPhi]pert[LI[n_]]->0}
translatedPert=perturbation/. {CD[a_]:>PD[a],H[LI[order_],a_,b_]:>If[order==1,r[a,b],0]}/. {s_Symbol/;StringMatchQ[SymbolName[s],"\[Tau]*"]:>ToExpression["d"<>StringDrop[SymbolName[s],1]],\[Alpha]->a,\[Beta]->b,\[Rho]->c,\[Tau]->d,\[CapitalPhi]->\[CurlyPhi]}//CollectTensors//ScreenDollarIndices
LagUV=FreeLag[\[Phi]]+(XActToMatchete[translatedPert,{r->h,\[CurlyPhi]->\[Phi]}]/.{\[Kappa]->\[Kappa][],M->M[]});
eftOrder = 7
Print["Matching, order "<>ToString[eftOrder]]
$printProgress=True
CheckAbort[LEFT0 = Match[LagUV, LoopOrder -> 1, EFTOrder -> eftOrder];
Print["Saving result"]
DumpSave["LEFT.mx", {LEFT0}],
Print["Abort detected"]]
