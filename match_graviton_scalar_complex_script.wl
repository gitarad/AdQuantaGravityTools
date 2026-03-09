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
DefineField[\[Phi], Matchete`Scalar, SelfConjugate -> False, 
  Mass -> {Heavy, M}];
DefineField[h, Graviton, SelfConjugate -> True, Mass -> {Light, 0}];
DefineCoupling[\[Kappa], SelfConjugate -> True];
Print[GetFields[]];
\[Phi]bar[]:=Bar[\[Phi][]];
DefManifold[Mpert,4,{\[Alpha],\[Beta],\[Rho],\[Tau]}];
DefMetric[{1,3,0},g[-\[Alpha],-\[Beta]],CD,{";","\[Del]"}];
DefMetricPerturbation[g,H,\[CurlyEpsilon]];
DefTensor[\[CapitalPhi][],Mpert];
DefTensor[\[CapitalPhi]c[],Mpert];
DefTensorPerturbation[\[CapitalPhi]pert[LI[0]],\[CapitalPhi][],\[CurlyEpsilon]];
DefTensorPerturbation[\[CapitalPhi]pertc[LI[0]],\[CapitalPhi]c[],\[CurlyEpsilon]];
DefManifold[M4,4,{a,b,c,d}];
DefMetric[{1,3,0},\[Eta][-a,-b],PD,{",","\[PartialD]"},FlatMetric->True];
DefTensor[r[-a,-b],M4,Symmetric[{-a,-b}]];
DefTensor[R[-a,-b],M4,Symmetric[{-a,-b}]];
DefTensor[\[CurlyPhi][],M4];
DefTensor[\[CurlyPhi]c[],M4];
DefConstantSymbol[{hbar,hbar2,\[Kappa],\[Kappa]2,\[Epsilon],M,\[Mu]bar2}];
Print["Expanding perturbation"]
perturbationOrder=3;
perturbationSum[expr_,order_,orderParam_]:=Sum[PerturbFlat[expr,i]orderParam^i/Factorial[i],{i,1,order}]
perturbation = (ExpandPerturbation[perturbationSum[Sqrt[-Detg[]](g[\[Alpha],\[Beta]]CD[-\[Alpha]][\[CapitalPhi][]]CD[-\[Beta]][\[CapitalPhi]c[]]-M^2\[CapitalPhi][]\[CapitalPhi]c[]),perturbationOrder,\[Kappa]]]//ToFlat//CollectTensors)/.{Sqrt[-Detg[]]->1,\[CapitalPhi]pert[LI[n_]]->0,\[CapitalPhi]pertc[LI[n_]]->0};
translatedPert=perturbation/. {CD[a_]:>PD[a],H[LI[order_],a_,b_]:>If[order==1,r[a,b],0]}/. {s_Symbol/;StringMatchQ[SymbolName[s],"\[Tau]*"]:>ToExpression["d"<>StringDrop[SymbolName[s],1]],\[Alpha]->a,\[Beta]->b,\[Rho]->c,\[Tau]->d,\[CapitalPhi]->\[CurlyPhi],\[CapitalPhi]c->\[CurlyPhi]c}//CollectTensors//ScreenDollarIndices;
LagUV=FreeLag[\[Phi]]+(XActToMatchete[translatedPert,{r->h,\[CurlyPhi]->\[Phi],\[CurlyPhi]c->\[Phi]bar}]/.{\[Kappa]->\[Kappa][],M->M[]});
Print["Matching"]
$dontCheckLagrangian=True;
LEFT0 = Match[LagUV, LoopOrder -> 1, EFTOrder -> 6];
Print["Saving result"]
DumpSave["LEFT.mx", {LEFT0}];
