(* Verification of the alternating Hurwitz zeta function.
   Run this file from the verification directory using Get
   or WolframKernel -script. The reference CSV is read from Dados.
   The Wolfram kernel performs the evaluations. *)

(* Alternating Hurwitz zeta function *)
(* The calculation uses the difference of two Hurwitz zeta functions. At s = 1, the digamma form avoids subtracting the two poles. *)
inicio = Now;
versao = $Version;
pasta = If[StringQ[$InputFileName] && StringEndsQ[$InputFileName,"verification.wl"],
   DirectoryName[ExpandFileName[$InputFileName]], NotebookDirectory[]];
projeto = pasta;
{versao, DateString[inicio,{"ISODate"," ","Time"," ","TimeZoneName"}]};

zetaAlternada[s_,q_] := If[TrueQ[s == 1],
   (PolyGamma[0,(q+1)/2]-PolyGamma[0,q/2])/2,
   2^(-s) (Zeta[s,q/2]-Zeta[s,(q+1)/2])];
FullSimplify[zetaAlternada[1,1] == Log[2]];

(* Normalization and transition regime *)
(* The parameter q is real, with q ≥ 100 and −1 ≤ x ≤ 1. The calculation reconstructs t from the rational coordinates at each working precision. *)
altura[q_,x_] := Pi q + x Sqrt[q];
valorZ[q_,x_,p_Integer] := Module[{t=N[altura[q,x],p]},
   Exp[I t Log[q]] zetaAlternada[1/2+I t,q]];

(* Functions F and G *)
(* The phase and the argument of the complementary error function follow the manuscript's conventions. *)
perfil[x_] := (1+I)/2 Exp[-I x^2/(2 Pi)] Erfc[-(1-I)x/(2 Sqrt[Pi])];
correcao[x_,f_] := 1/2+I/(6 Pi)-x^2 (1-I x f)/(6 Pi^2);
{perfil[0],correcao[0,perfil[0]]};

(* Two-term approximation *)
(* The first term is F. The correction added to it is G/√q. *)
doisTermos[q_,f_,g_] := f + g/Sqrt[q];

(* Absolute and scaled errors *)
(* The calculation records the real and imaginary parts, the one-term error, |R|, and q|R|. The working precision applies to the evaluation. Subtraction may reduce the precision of the remainder. *)
calcular[q_,x_,p_Integer] := Module[{z,f,g,r},
   z=valorZ[q,x,p]; f=N[perfil[x],p]; g=N[correcao[x,f],p];
   r=z-doisTermos[q,f,g];
   {q,x,p,Re[z],Im[z],Re[f],Im[f],Re[g],Im[g],Re[r],Im[r],Abs[r],q Abs[r],Abs[z-f]}];
campos={"q","x","working_dps","Z_re","Z_im","F_re","F_im","G_re","G_im","R_re","R_im","abs_R","q_abs_R","abs_one_term"};

(* Symbolic identities *)
(* Each line simplifies a difference to zero. The assumptions for contour rotation and integration of the remainder are part of the analytical proof. *)
simbolicos=Module[{a,v,w,y,h,x,z,r,k,theta,delta,expoente,fase,p1,thetaY,difase,ff,f1,f2,f3,derivar,raw},
   k=(a w-v)/(1-a); theta=ArcTan[r/(Sqrt[2]+r)]; delta=r/Sqrt[2]-theta;
   expoente=-Pi z/h-(1/2+I Pi/h^2+I x/h) Log[1+I h z];
   fase=x z-I Pi z^2/2; p1=-I z/2-I x z^2/2-Pi z^3/3;
   thetaY=ArcTan[y h^2]; difase=2 Pi y-2 (Pi/h^2+x/h)thetaY+I thetaY;
   f1=(1-I x ff)/Pi; derivar[expr_]:=D[expr,x]+f1 D[expr,ff];
   f2=derivar[f1]; f3=derivar[f2];
   raw=<|"K and the splitting of the integral"->((w-v)/(1-a)-w-k),
   "Exponential form of K"->(k-w (1-v/(a w))/(1/a-1)),
   "K at the origin"->(- (D[difase,y]/.y->0)/(2 Pi)-x h/Pi+I h^2/(2 Pi)),
   "Taylor expansion of the phase, terms 0 and 1"->(Normal[Series[expoente,{h,0,1}]]-fase-h p1),
   "Derivative of the rotated argument"->(D[theta,r]-1/(Sqrt[2](1+Sqrt[2]r+r^2))),
   "Derivative of delta'/r"->(D[D[delta,r]/r,r]+(1/Sqrt[2]+2r+r^2/Sqrt[2])/(1+Sqrt[2]r+r^2)^2),
   "Differential equation for F"->(Pi D[perfil[x],x]-1+I x perfil[x]),
   "F at the origin"->(perfil[0]-(1+I)/2),
   "Reduction of the coefficient G"->(-I f1/2-I x f2/2-Pi f3/3-(correcao[x,ff]-1/2))|>;
   FullSimplify[#,Assumptions->(r>0 && 0<a<1 && w!=0 && Element[x,Reals])]& /@raw];
Grid[Prepend[KeyValueMap[List,simbolicos],{"Identity","Difference"}],Alignment->Left];

(* Upper bounds and rational arithmetic *)
(* The recurrence for even moments starts at 67/80, an upper bound on the zeroth Gaussian moment. The sum is computed from the polynomial coefficients without rounding. *)
piInf=157/50; piSup=22/7; hMax=1/10; aGauss=9/8; kExp=67/20;
momentoParidade[j_] := FullSimplify[j! Zeta[j+1]/(2 Pi)^(j+1)];
momentoImpar[j_] := FullSimplify[j!(1-2^(-j-1)) Zeta[j+1]/Pi^(j+1)];
cotaPartida[m_] := m ((2+hMax) momentoParidade[1]+2 piSup hMax^3 momentoParidade[3]/3)+1/100000;
bS=cotaPartida[5]; bTres=cotaPartida[3];
bGlobal=(2+hMax) momentoImpar[1]+2 piSup hMax^3 momentoImpar[3]/3;
polinomio=Expand[(u/2+u^2/2+piSup u^3/3)^2+u^2/2+2u^3/3+piSup u^4/2];
gauss[0]=67/80; gauss[1]=1/(2 aGauss);
gauss[j_Integer /; j>=2] := (j-1)gauss[j-2]/(2 aGauss);
jBarra=Total[Table[Coefficient[polinomio,u,j](7 gauss[j]/4+j!/kExp^(j+1)),{j,2,6}]];
{polinomio,bS,bGlobal,bTres,jBarra};

exInf[v_,n_] := Sum[v^k/k!,{k,0,n}];
exSup[v_,n_] := exInf[v,n]+v^(n+1)/((n+1)! (1-v/(n+2)));
atanFaixa[v_,n_] := With[{s=Sum[(-1)^k v^(2k+1)/(2k+1),{k,0,n-1}]},
   Sort[{s,s+(-1)^n v^(2n+1)/(2n+1)}]];
r0=3/10; raizInf=140/99; raizSup=99/70;
v0=r0/(r0+raizInf); d0=r0/raizSup-v0+v0^3/3-v0^5/5;
passo=r0/4096;
dSoma=passo Total[Table[With[{z=(1414213/10^6) j passo+(j passo)^2},
   Floor[10^12 z/((1414214/10^6)(1+z))]],{j,0,4095}]]/10^12;
transicao[d_] := piInf d/r0^2>9/8+5/22 && 10 piInf d/r0-1/raizInf>67/20;
argMax=(2hMax+hMax^2)(461/100)+2 piSup(461/100)^3 hMax^4/3;
contasExatas=<|"Bounds on pi (Machin)"->With[{v=atanFaixa[1/5,6],w=atanFaixa[1/239,4]},piInf<16 First[v]-4 Last[w] && 16 Last[v]-4 First[w]<piSup],
"Square-root bounds and transition"->(raizInf^2<2<raizSup^2 && (1414213/10^6)^2<2<(1414214/10^6)^2 && transicao[d0] && transicao[dSoma] && dSoma==1592131749274149/40960000000000000),
"Kernel moments"->({momentoParidade[1],momentoParidade[3],momentoImpar[1],momentoImpar[3]}=={1/24,1/240,1/8,1/16}),
"Bounding polynomial"->(Rest[Rest[CoefficientList[polinomio,u]]]=={3/4,7/6,241/84,22/21,484/441}),
"Gaussian and exponential bounds"->(piSup/(4aGauss)<(67/80)^2 && (1/2)/(4(5/22))==11/20 && exSup[11/20,6]<7/4),
"Cutoff and local multiplier"->(exInf[461/100,16]>100 && exSup[2,12]<100 && exSup[461/1000,6]<5/3 && 5/3(exSup[1,6]-1)<3 && argMax==51902943991/52500000000 && argMax<989/1000),
"Tail"->(piInf>3 && (1/100000)^2==(1/100)^5 && 100^(-6)<1/3),
"B_s, B_g and J"->({bS,bGlobal,jBarra}=={689147/1575000,22061/84000,533948559110614545113/69278297361197722560})|>;
Grid[KeyValueMap[List,contasExatas],Alignment->Left];

(* Comparison of 4.3/q and 4.117/q *)
(* The main result uses a multiplier of 5 on the short interval. The appendix uses 3. The global kernel estimate gives a third calculation. *)
cPrincipal=jBarra/2+bS; cGlobal=jBarra/2+bGlobal; cApendice=jBarra/2+bTres;
somas=<|"Main bound"->cPrincipal,"Global kernel bound"->cGlobal,"Appendix bound"->cApendice|>;
estritas=<|"Main bound as a fraction"->(cPrincipal==371609064358084350147257/86597871701497153200000),
"Global bound as a fraction"->(cGlobal==14258445238530472531037/3463914868059886128000),
"Appendix bound as a fraction"->(cApendice==356452924783991052355217/86597871701497153200000),
"Margin below 43/10"->(43/10-cPrincipal==761783958353408612743/86597871701497153200000 && cPrincipal<43/10),
"Both alternatives below 4117/1000"->(cApendice<cGlobal<4117/1000)|>;
Grid[Prepend[KeyValueMap[{#1,#2,N[#2,15]}&,somas],{"Calculation","Exact value","Decimal value"}],Alignment->Left];

(* Numerical tests and comparison with the original results *)
(* Ten quantities are compared in each row, using the coordinates and working precision as the key. The CSV values serve only as a reference. The manuscript table is rounded to twelve decimal places. *)
qs={100,100001/1000,201/2,101,150,400,1600};
xs=Range[-1,1,1/20];
linhas55=Flatten[Table[calcular[q,x,55],{q,qs},{x,xs}],1];
porQ=GatherBy[linhas55,First];
maximos=First[MaximalBy[#,Function[linha,linha[[13]]]]]& /@ porQ;
linhas75=calcular[#[[1]],#[[2]],75]& /@ maximos;
todas=Join[linhas55,linhas75];
{Length[linhas55],Length[linhas75]};

decimalExato[str_String] := Module[{s=ToLowerCase[StringTrim[str]],partes,mant,m,casas,exp=0},
   If[!StringMatchQ[s,RegularExpression["[+-]?[0-9]+(?:\\.[0-9]*)?(?:e[+-]?[0-9]+)?"]],Return[$Failed]];
   partes=StringSplit[s,"e"]; mant=First[partes]; If[Length[partes]==2,exp=ToExpression[Last[partes]]];
   m=StringTrim[mant,StartOfString~~("+"|"-")]; casas=If[StringContainsQ[m,"."],StringLength[m]-First[First[StringPosition[m,"."]]],0];
   If[StringStartsQ[mant,"-"],-1,1] FromDigits[StringDelete[m,"."]] 10^(exp-casas)];
diferenca[a_,b_] := Abs[SetPrecision[a,Infinity]-SetPrecision[b,Infinity]];
exigir[cond_,msg_String] := If[!TrueQ[cond],Print[msg];Abort[]];
csvOrig=FileNameJoin[{projeto,"Dados","Hurwitz_Numerical_Results.csv"}];
texOrig=FileNameJoin[{projeto,"Manuscrito","HurwitzZetaManuscript.tex"}];
nbOrig=FileNameJoin[{projeto,"Wolfram","Hurwitz_Calculations.nb"}];
exigir[And@@(FileExistsQ/@{csvOrig,texOrig,nbOrig}),"A reference file is missing from the project folder."];
csv=Import[csvOrig,"CSV","Numeric"->False];
exigir[Take[First[csv],14]==Join[{"set"},Take[campos,13]],"The CSV header differs from the expected header."];
orig=Select[Rest[csv],First[#]=="manuscript" && MemberQ[{"55","75"},#[[4]]]&];
ref=Association[Map[(decimalExato /@ #[[2;;4]]) -> (decimalExato /@ #[[5;;14]])&,orig]];
chaves=#[[1;;3]]&/@todas;
exigir[Length[orig]==294 && Sort[Keys[ref]]===Sort[chaves] && FreeQ[ref,$Failed],"The CSV contains invalid (q, x, precision) keys or decimal values."];
desvios=Map[Function[linha,MapThread[diferenca,{linha[[4;;13]],ref[linha[[1;;3]]]}]],todas];

texto=Import[texOrig,"Text"];
parteTabela=First[StringSplit[Last[StringSplit[texto,"\\label{tab:max}"]],"\\end{tabular}"]];
tabImpressa=(decimalExato/@#)&/@StringCases[parteTabela,RegularExpression["(?m)^\\s*([0-9.]+)\\s*&\\s*(-?[0-9.]+)\\s*&\\s*([0-9.]+)\\s*\\\\\\\\"] -> {"$1","$2","$3"}];
exigir[Dimensions[tabImpressa]=={7,3} && tabImpressa[[All,1]]===qs,"The seven rows of tab:max were not found in the expected order."];
nbLido=ToExpression[Import[nbOrig,"Text"],InputForm,HoldComplete];
curvas=Cases[nbLido,LineBox[pts_List] /; Dimensions[pts]=={41,2} && Max[Abs[pts[[All,1]]-xs]]<10^-14 :> pts,Infinity];
exigir[Length[curvas]==7,"The plot in the original notebook does not contain the seven expected curves."];
difCurvas=MapThread[diferenca,{linhas55[[All,13]],Flatten[curvas[[All,All,2]]]}];
difImpresso=MapThread[diferenca,{maximos[[All,13]],tabImpressa[[All,3]]}];
difPrecisao=MapThread[diferenca,{maximos[[All,13]],linhas75[[All,13]]}];
rotulos=campos[[4;;13]];
maxPorCampo=Transpose[{rotulos,Max/@Transpose[Take[desvios,287]],Max/@Transpose[Take[desvios,-7]]}];
Grid[Prepend[({#[[1]],N[#[[2]],5],N[#[[3]],5]}& /@maxPorCampo),{"Field","Maximum difference, 55 digits","Maximum difference, 75 digits"}],Alignment->Left];

quadratura[q_,x_] := Module[{s=1/2+I altura[q,x],f,y},
   f[v_?NumericQ]:=If[TrueQ[v==0],-I s/(Pi q),((1+I v/q)^(-s)-(1-I v/q)^(-s))/(2 Sinh[Pi v])];
   1/(2 Sqrt[q])+I/Sqrt[q] NIntegrate[f[y],{y,0,1,Sqrt[q],4 Sqrt[q],Infinity},
   WorkingPrecision->60,AccuracyGoal->35,PrecisionGoal->35,MaxRecursion->30]];
difAP=Table[diferenca[quadratura[100,x],valorZ[100,x,70]],{x,{-1,0,1}}];
testes=<|"s=1, q=1"->FullSimplify[zetaAlternada[1,1]==Log[2]],"Nine identities"->(Values[simbolicos]===ConstantArray[0,9]),
"Exact upper bounds"->And@@Values[contasExatas],"Sums and margins"->And@@Values[estritas],
"287 + 7 keys"->(Length[linhas55]==287 && Length[linhas75]==7 && Length[DeleteDuplicates[chaves]]==294),
"CSV, 55 digits (10 fields)"->(Max[Take[desvios,287]]<10^-44),
"CSV, 75 digits (10 fields)"->(Max[Take[desvios,-7]]<10^-64),
"Stability of q|R| at 55 and 75 digits"->(Max[difPrecisao]<10^-44),
"Manuscript table (12 decimal places and x)"->(Max[difImpresso]<=1/(2 10^12) && maximos[[All,2]]===tabImpressa[[All,2]]),
"Original curves (machine precision)"->(Max[difCurvas]<10^-15),
"Abel–Plana at three points"->(Max[difAP]<10^-35),
"Main bound at the sampled points"->(Max[todas[[All,13]]]<43/10),
"Appendix bound at the sampled points"->(Max[todas[[All,13]]]<4117/1000)|>;
Grid[KeyValueMap[{#1,If[TrueQ[#2],"Passed","Discrepancy"]}&,testes],Alignment->Left];

(* Two-dimensional plots *)
(* Line segments connect the grid points. The continuous bound curves are computed from the indicated formulas. *)
estilo={FontFamily->"Times",FontSize->12,FontColor->Black};
estilosCurvas={Directive[Black,Thick],Directive[GrayLevel[0.35],Dashed,Thick]};
grafMax=ListLinePlot[Transpose[{N[qs],N[maximos[[All,13]]]}],Frame->True,Axes->False,
   FrameLabel->{"q","Observed maximum of q |R|"},ScalingFunctions->{"Log",None},
   PlotStyle->Black,PlotMarkers->Automatic,PlotRange->All,BaseStyle->estilo,ImageSize->480];
grafMax;

corte=First[porQ];
grafTermos=ListLinePlot[{corte[[All,{2,14}]],corte[[All,{2,12}]]},ScalingFunctions->{None,"Log"},
   PlotStyle->estilosCurvas,PlotLegends->{"F only","F + G/√q"},Frame->True,Axes->False,
   FrameLabel->{"x (q = 100)","Absolute error"},PlotRange->All,BaseStyle->estilo,ImageSize->470];
grafTermos;

grafCotas=Plot[{43/(10 q),4117/(1000 q)},{q,100,1600},ScalingFunctions->{"Log","Log"},
   PlotStyle->estilosCurvas,PlotLegends->{"4.3/q","4.117/q"},Frame->True,Axes->False,
   FrameLabel->{"q","Bound"},BaseStyle->estilo,ImageSize->470];
grafCotas;

(* Three-dimensional plots *)
(* The first surface connects the 287 values of |R|. The second interpolates q|R| linearly in x and log₁₀(q). Heights between the grid points are interpolated. *)
pontos3D=({#[[2]],Log[10,#[[1]]],#[[12]]}& /@linhas55);
grafSuperficie=ListPlot3D[N[pontos3D],Mesh->All,MeshStyle->Directive[GrayLevel[0.4],Thin],
   PlotStyle->GrayLevel[0.8],Lighting->"Neutral",AxesLabel->{"x","log₁₀ q","|R|"},
   PlotRange->All,BoxRatios->{1.4,1,0.65},BaseStyle->estilo,ImageSize->480];
grafSuperficie;

interp=ListInterpolation[N[porQ[[All,All,13]]],{N[Log[10,qs]],N[xs]},InterpolationOrder->1];
grafEscalado=Plot3D[interp[logq,x],{x,-1,1},{logq,2,N[Log[10,1600]]},
   PlotPoints->{41,17},MaxRecursion->0,Mesh->{8,6},MeshStyle->GrayLevel[0.4],PlotStyle->GrayLevel[0.75],
   Lighting->"Neutral",AxesLabel->{"x","log₁₀ q","q |R|"},PlotRange->All,
   BoxRatios->{1.4,1,0.65},BaseStyle->estilo,ImageSize->480];
grafEscalado;

(* Complete results tables *)
(* Each of the seven tables contains all 41 values of x. The CSV also records the complex components and the results of the seven runs repeated at 75-digit precision. The comparison uses absolute tolerances. *)
cabMax={"q","x","q |R| at 55 digits","q |R| at 75 digits","Difference between 55 and 75 digits"};
Grid[Prepend[MapThread[{#1[[1]],#1[[2]],N[#1[[13]],16],N[#2[[13]],16],N[#3,5]}&,
  {maximos,linhas75,difPrecisao}],cabMax],Alignment->Left];

tabelaSerie[filas_] := Grid[Prepend[({#[[2]],N[#[[14]],11],N[#[[12]],11],N[#[[13]],12]}& /@filas),
   {"x","|Z − F|","|R|","q |R|"}],Alignment->Right,Dividers->{None,{False,True}},Spacings->{1.4,0.45}];

tabelaSerie[porQ[[1]]];

tabelaSerie[porQ[[2]]];

tabelaSerie[porQ[[3]]];

tabelaSerie[porQ[[4]]];

tabelaSerie[porQ[[5]]];

tabelaSerie[porQ[[6]]];

tabelaSerie[porQ[[7]]];

numeroCSV[v_?NumberQ] := Module[{nd,d,e,seq,m},
   If[TrueQ[v==0],Return["0"]]; If[IntegerQ[v],Return[ToString[v]]];
   nd=If[Precision[v]===Infinity,16,Max[1,Floor[Precision[v]]]];
   {d,e}=RealDigits[Abs[SetPrecision[v,Infinity]],10,nd];
   seq=StringJoin[ToString/@d]; m=StringTake[seq,1]<>If[nd>1,"."<>StringDrop[seq,1],""];
   If[v<0,"-",""]<>m<>If[e==1,"","e"<>ToString[e-1]]];
outCsv=FileNameJoin[{pasta,"verification_results.csv"}];
csvDados=MapThread[Join[numeroCSV/@#1,numeroCSV/@{Max[#2],Last[#2]}]&,{todas,desvios}];
Export[outCsv,Prepend[csvDados,Join[campos,{"max_abs_diff_reference_10_fields","abs_diff_q_abs_R_reference"}]],"CSV"];
fim=Now;
estado=If[And@@Values[testes],"All tests passed","Discrepancies found. See the test results."];
{estado,Length[todas],DateString[fim,{"ISODate"," ","Time"," ","TimeZoneName"}]};
Print[estado];
