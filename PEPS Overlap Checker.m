(* ::Package:: *)

(* ::Input:: *)
(*Eprod[a_]:=Flatten[a.Conjugate[Transpose[a,{2,3,4,5,1}]],{{1,5},{2,6},{3,7},{4,8}}]*)
(*Eprod[a_,b_]:=Flatten[a.Conjugate[Transpose[b,{2,3,4,5,1}]],{{1,5},{2,6},{3,7},{4,8}}]*)


(* ::Input:: *)
(*Lprod[L_,R_]:=Flatten[Transpose[L,{1,4,2,3}].R,{{1},{4},{2,5},{3,6}}]*)
(*Vprod[U_,D_]:=Flatten[U.Transpose[D,{2,3,1,4}],{{1,4},{2,5},{3},{6}}]*)


(* ::Input:: *)
(*line[Lup_,Ldo_]:=Module[{acc,er},*)
(*acc=Eprod[Lup[[1]],Ldo[[1]]];*)
(*Do[*)
(*er=Eprod[Lup[[i+1]],Ldo[[i+1]]];*)
(*acc=Lprod[acc,er];*)
(*,{i,1,3}];*)
(*acc]*)
(**)


(* ::Input:: *)
(*Contract[Pup_,Pdo_]:=Module[{myline,acc},*)
(*acc=line[Pup[[All,4]],Pdo[[All,4]]];*)
(*Do[*)
(*myline=line[Pup[[All,i-1]],Pdo[[All,i-1]]];*)
(*acc=Vprod[acc,myline];*)
(*,{i,4,2,-1}];*)
(*Flatten[acc][[1]]*)
(*]*)


(* ::Input:: *)
(*func[l_,r_,u_,d_,s_,x_,y_]:= I*(s-1)+l*3.0+(r/(x^2+I Sqrt[y]))^2+Sqrt[u]+I*Log[d](* I(s-1)+l+r+u+d-4 *)*)
(*SymTensor[L_,R_,U_,D_,x_,y_]:=Module[{A,a,S=2},*)
(*A=Array[a,{L,R,U,D,2}];*)
(*Do[*)
(*Do[*)
(*Do[*)
(*Do[*)
(*Do[*)
(*A[[l,r,u,d,s]]=func[l,r,u,d,s,x,y]*)
(*,{s,1,S}]*)
(*,{d,1,D}]*)
(*,{u,1,U}]*)
(*,{r,1,R}]*)
(*,{l,1,L}];*)
(*A*)
(*]*)


(* ::Input:: *)
(*Peps1=Array[f,{4,4}];*)
(*Peps2=Array[f,{4,4}];*)
(*Do[*)
(*Do[*)
(*Peps1[[x,y]]=SymTensor[2,2,2,2,x,y];*)
(*Peps2[[x,y]]=(1/2) SymTensor[2,2,2,2,x,y];*)
(*,{x,2,3}]*)
(*,{y,2,3}];*)
(*Module[{x,y},*)
(*Do[*)
(*x=1;y=n;*)
(*Peps1[[x,y]]=SymTensor[1,2,2,2,x,y];*)
(*Peps2[[x,y]]=(1/4) SymTensor[1,2,2,2,x,y];*)
(*x=4;y=n;*)
(*Peps1[[x,y]]=SymTensor[2,1,2,2,x,y];*)
(*Peps2[[x,y]]=(1/4) SymTensor[2,1,2,2,x,y];*)
(*x=n;y=1;*)
(*Peps1[[x,y]]=SymTensor[2,2,2,1,x,y];*)
(*Peps2[[x,y]]=(1/4) SymTensor[2,2,2,1,x,y];*)
(*x=n;y=4;*)
(*Peps1[[x,y]]=SymTensor[2,2,1,2,x,y];*)
(*Peps2[[x,y]]=(1/4) SymTensor[2,2,1,2,x,y];*)
(**)
(*,{n,2,3}];*)
(*x=1;y=1;*)
(*Peps1[[x,y]]=SymTensor[1,2,2,1,x,y];*)
(*Peps2[[x,y]]=3 SymTensor[1,2,2,1,x,y];*)
(*x=4;y=1;*)
(*Peps1[[x,y]]=SymTensor[2,1,2,1,x,y];*)
(*Peps2[[x,y]]=3 SymTensor[2,1,2,1,x,y];*)
(*x=1;y=4;*)
(*Peps1[[x,y]]=SymTensor[1,2,1,2,x,y];*)
(*Peps2[[x,y]]=3 SymTensor[1,2,1,2,x,y];*)
(*x=4;y=4;*)
(*Peps1[[x,y]]=SymTensor[2,1,1,2,x,y];*)
(*Peps2[[x,y]]=3 SymTensor[2,1,1,2,x,y];*)
(*]*)
(**)


(* ::Input:: *)
(*N[Contract[Peps1,Peps2],20]*)


(* ::Input:: *)
(*1.591182839737796`*^38- 9.271418372179549`*^20 I*)


(* ::Input:: *)
(*Abs[%]*)


(* ::Input:: *)
(*norm=0.612943975137486;*)
(*can=0.661086186100242;*)


(* ::Input:: *)
(*Sqrt[norm]*)
(*Sqrt[can]*)


(* ::Input:: *)
(*normal=3.59295799358889*)


(* ::Input:: *)
(*direct=3.56911300722999/Sqrt[3.16460664004381 9.32102735763028]*)


(* ::Input:: *)
(*normal/direct*)


(* ::Input:: *)
(*direct*Sqrt[3.16460664004381 9.32102735763028]*)


(* ::InheritFromParent:: *)
(*3.56911300722999`*)


(* ::Input:: *)
(*Sqrt[10.0747051431802]*)


(* ::Input:: *)
(*3.61731382401674/Sqrt[3.17406760217551]*)


(* ::Input:: *)
(*%/Sqrt[9.32102735763027]*)


(* ::Input:: *)
(* 3.61731382401674 *)


(* ::Input:: *)
(*0.361170763674116*9.32102735763027*)


(* ::Input:: *)
(*0.665037615980690*(9.32102735763027*3.17406760217551)*)
