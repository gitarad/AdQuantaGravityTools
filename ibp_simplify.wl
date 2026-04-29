(* ::Package:: *)

Print["Loading xAct"]
Needs["xAct`xTensor`"];
Needs["xAct`xPert`"];
Needs["xAct`xTras`"];
DefManifold[M4,4,{a,b,c,d}]
DefMetric[{1,3,0},\[Eta][-a,-b],PD,{",","\[PartialD]"},FlatMetric->True]
DefTensor[r[-a,-b],M4,Symmetric[{-a,-b}]]
DefConstantSymbol[{hbar,hbar2,\[Kappa],\[Kappa]2,\[Epsilon],M,\[Mu]bar2,G}]
DefManifold[Mpert,4,{\[Alpha],\[Beta],\[Rho],\[Tau]}]
DefMetric[{1,3,0},g[-\[Alpha],-\[Beta]],CD,{";","\[Del]"}]
DefMetricPerturbation[g,H,\[CurlyEpsilon]]
DefTensor[R[-a,-b],M4,Symmetric[{-a,-b}]]
DefConstantSymbol[{A1,A2,A3}]
IBPSimplify[expr_,varSym_,dumpName_]:=Module[
{beforeAnsatz,ansatz,EquationOfConstraint,ReducedEq,eqs,rhsConstants,rhsSubstitutions,ibpZeroOperator},
Print["Collecting tensors"];
beforeAnsatz = List@@(expr//CollectTensors);
ansatz=MakeAnsatz[beforeAnsatz];
Export["IBPSimplifyAnsatz"<>dumpName<>".bin",BinarySerialize[ansatz],"binary"];
Print["Calculating equation"];
EquationOfConstraint=VarD[varSym][ansatz];
Print["Simplifying equation"];
ReducedEq=CollectTensors[EquationOfConstraint//ScreenDollarIndices//ToCanonical//ContractMetric//SortCovDs//Simplify];
Print["Dumping equation"];
Export["IBPSimplifyEqs"<>dumpName<>".bin",BinarySerialize[ReducedEq],"binary"]
Print["Solving equation"];
eqs=SolveConstants[ReducedEq==0,!{hbar,hbar2,\[Kappa],\[Kappa]2,\[Epsilon],M,\[Mu]bar2,G,A1,A2,A3}][[1]];
Print["substituting solution"];
rhsConstants=Flatten[Table[Cases[eq[[2]],s_Symbol/;StringMatchQ[SymbolName[s],"C"~~DigitCharacter..],Infinity]//DeleteDuplicates,{eq,eqs}]]//DeleteDuplicates;
rhsSubstitutions=Table[Cnum->-1,{Cnum,rhsConstants}];
ibpZeroOperator=ansatz/.eqs/.rhsSubstitutions;
ibpZeroOperator=ibpZeroOperator/.s_Symbol/;StringMatchQ[SymbolName[s],"C"~~DigitCharacter..]:>-1;
(expr+ibpZeroOperator)//CollectTensors
];
param = $CommandLine[[4]];
Print["Parameter: ",param];
Print["Loading file"];
Get["IBPSimplifyInput"<>param<>".mx"];
Print["Simplifying"];
IBPSimplifyOutput=IBPSimplify[IBPSimplifyInput,r[-a,-b],param];
Print["Collecting tensors"];
IBPSimplifyOutput=IBPSimplifyOutput//CollectTensors;
Print["Saving result"];
DumpSave["IBPSimplifyOutput"<>param<>".mx", {IBPSimplifyOutput}];
Print["Saved result"];
