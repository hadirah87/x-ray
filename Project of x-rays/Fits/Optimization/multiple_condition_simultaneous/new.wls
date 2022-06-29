#!/usr/bin/env wolframscript
(* ::Package:: *)

(*Plot[UnitStep[x+45]-UnitStep[x-45],{x,-70,70}]*)


(*variables should be determined here*)
(*d=139   (*grating periodic distance of silica*)*)
variables=ToExpression[Rest[$ScriptCommandLine]]
d=Part[variables,1]
rho=Part[variables,6]
Subscript[\[Rho], sil]=rho*10^-4 (*scattering length density of silica*)
Subscript[\[Rho], coat]=1.0*Subscript[\[Rho], sil]
g=2*\[Pi]/d (*resiprocal lattice length*)
fd=Part[variables,2]
f=fd/100(*filling factor of grating*)
fbd=Part[variables,3]
fb=fbd/100
ftd=Part[variables,4]
ft=ftd/100
t=Part[variables,5]  (*tickness of grating *)
tc=t/2(*tickness of coated layer*)
Subscript[f, new]=0.0 (*filling factor of coated layer*)


j=15;
<<FourierSeries`
Quiet[Table[Subscript[\[Rho]c, m]=Re[NFourierCoefficient[Subscript[\[Rho], sil]*UnitStep[Cos[g*y/(2*fb)]],y,m, FourierParameters->{1, g} ]],{m,-2j,2j}]];
(*Subscript[\[Rho]c, 4]*)
(*\[Rho]1[y_,j_] := Sum[Subscript[\[Rho]c, m]*Cos[g*y*m],{m,-j,j}]*)



j=15;
<<FourierSeries`
(*Quiet[Table[Subscript[\[Rho]s, m]=Re[NFourierCoefficient[Subscript[\[Rho], coat]*UnitStep[Cos[g*y/(2*f)]],y,m, FourierParameters\[Rule]{1, g} ]],{m,-2j,2j}]];*)
Quiet[Table[Subscript[\[Rho]s, m]=Re[NFourierCoefficient[Subscript[\[Rho], sil]*UnitStep[Cos[g*y/(2*f)]],y,m, FourierParameters->{1, g} ]],{m,-2j,2j}]];
(*Subscript[\[Rho]s, 3]
\[Rho]3[y_,j_] := Sum[Subscript[\[Rho]s, m]*Cos[g*y*m],{m,-j,j}]*)



j=15;
<<FourierSeries`
Quiet[Table[Subscript[\[Rho]b, m]=Re[NFourierCoefficient[Subscript[\[Rho], sil]*UnitStep[Cos[g*y/(2*ft)]],y,m, FourierParameters->{1, g} ]],{m,-2j,2j}]];
(*Subscript[\[Rho]b, 4]
\[Rho]5[y_,j_] := Sum[Subscript[\[Rho]b, m]*Cos[g*y*m],{m,-j,j}]*)



(*j=9;
<<FourierSeries`
Quiet[Table[Subscript[\[Rho], n]=Re[NFourierCoefficient[rho*10^-4*UnitStep[Cos[g*x/(2*f)]], x,n, FourierParameters -> {1,2*Pi/140}]],{n,-2j,2j}]];*)


(*note : In mathematica 6 use the following command*)


(*//j=5;
<<FourierSeries`
Quiet[Table[Subscript[\[Rho], n]=Re[NFourierCoefficient[2.06*10^-4*UnitStep[Cos[\[Pi]*x/70]], x,n, FourierParameters -> {-1,1/140}]],{n,-2j,2j}]];//*)


(*Subscript[\[Theta], 0]=0.2*\[Pi]/180;*)
(*d=140;*)
(*t=65;*)
(*f=fd/10;*)


j=15;
Subscript[\[Theta], 0]=0.20*\[Pi]/180;
For[h=1;\[Lambda]=0.077012,\[Lambda] <=0.077012,\[Lambda] =\[Lambda] +0.001,
For[s=1;\[Phi]=-10.0*\[Pi]/180,\[Phi]<=10.0*\[Pi]/180,\[Phi]=\[Phi]+0.0451*\[Pi]/180;
Subscript[k, y,s,h]=Subscript[k, iy]=(2\[Pi])/\[Lambda]*Sin[\[Phi]];
Subscript[k, z,s,h]=Subscript[k, iz]=(2\[Pi])/\[Lambda]*Sin[Subscript[\[Theta], 0]];
Subscript[newkz, s,h]=((2 \[Pi])/\[Lambda])*Sin[Subscript[\[Theta], 0]];
Subscript[phi, s]=\[Phi];
Subscript[newky, s,h]=(2\[Pi])/\[Lambda]*Sin[\[Phi]];
(*Schrodinger equation in the region I*)
F1[n_,m_]:=(-(Subscript[k, iy]+n*g)^2+Subscript[k, iz]^2+Subscript[k, iy]^2)*KroneckerDelta[n,m];
V1[n_,m_]:=-4*\[Pi]*Subscript[\[Rho]c, n-m];
M1=Table[F1[n,m]+V1[n,m],{n,-j,j},{m,-j,j}];
{vals1,vecs1}=Eigensystem[M1];


(*Schrodinger equation in region*)
F2[n_,m_]:=(-(Subscript[k, iy]+n*g)^2+Subscript[k, iz]^2+Subscript[k, iy]^2)*KroneckerDelta[n,m];
V2[n_,m_]:=-4*\[Pi]*Subscript[\[Rho]s, n-m];
M2=Table[F2[n,m]+V2[n,m],{n,-j,j},{m,-j,j}];
{vals2,vecs2}=Eigensystem[M2];


(*Schrodinger equation in region III*)
F3[n_,m_]:=(-(Subscript[k, iy]+n*g)^2+Subscript[k, iz]^2+Subscript[k, iy]^2)*KroneckerDelta[n,m];
V3[n_,m_]:=-4*\[Pi]*Subscript[\[Rho]b, n-m];
M3=Table[F3[n,m]+V3[n,m],{n,-j,j},{m,-j,j}];
{vals3,vecs3}=Eigensystem[M3];


(*Wave vector in different region based on conservation energy*)
Table[(*Subscript[p0, s]=*) Subscript[X, s,h,n]=Subscript[x, n]=Sqrt[Subscript[k, iy]^2+Subscript[k, iz]^2-(Subscript[k, iy]+n*g)^2],{n,-j,j}];
Table[(*Subscript[p1, s]=*)Subscript[kz1, s,h,n]=Subscript[k1, z,n]=Sqrt[vals1[[n+j+1]]],{n,-j,j}];
Table[(*Subscript[p2, s]=*)Subscript[kz2, s,h,n]=Subscript[k2, z,n]=Sqrt[vals2[[n+j+1]]],{n,-j,j}];
Table[(*Subscript[p3, s]=*)Subscript[kz3, s,h,n]=Subscript[k3, z,n]=Sqrt[vals3[[n+j+1]]],{n,-j,j}];
Table[Subscript[Pt, s,h,n]=Subscript[p, n]=Sqrt[Subscript[k, iy]^2+Subscript[k, iz]^2-4*\[Pi]*Subscript[\[Rho], sil]-(Subscript[k, iy]+n*g)^2],{n,-j,j}];
(*And the eigenvector coefficients are:*)
Table[Subscript[b1, n,m]=vecs1[[n+j+1,m+j+1]],{n,-j,j},{m,-j,j}];
Table[Subscript[b2, n,m]=vecs2[[n+j+1,m+j+1]],{n,-j,j},{m,-j,j}];
Table[Subscript[b3, n,m]=vecs3[[n+j+1,m+j+1]],{n,-j,j},{m,-j,j}];

(*Applying series of boundary condition*)
Arr={Table[KroneckerDelta[m,0]+Subscript[r, 0,m]-Sum[(Subscript[t, 1,n]+Subscript[r, 1,n]*E^(I*Subscript[k1, z,n]*tc))*Subscript[b1, n,m],{n,-j,j}]==0,{m,-j,j}],Table[-Subscript[k, iz]*KroneckerDelta[m,0]+Subscript[x, m]*Subscript[r, 0,m]+Sum[Subscript[k1, z,n]*(Subscript[t, 1,n]-Subscript[r, 1,n]*E^(I*Subscript[k1, z,n]*tc))*Subscript[b1, n,m],{n,-j,j}]==0,{m,-j,j}],

Table[Sum[(Subscript[r, 1,n]+Subscript[t, 1,n]*E^(I*Subscript[k1, z,n]*tc))*Subscript[b1, n,m],{n,-j,j}]-Sum[(Subscript[t, 2,n]+Subscript[r, 2,n]*E^(I*Subscript[k2, z,n]*(t-tc)))*Subscript[b2, n,m],{n,-j,j}]==0,{m,-j,j}],Table[Sum[Subscript[k1, z,n]*(Subscript[r, 1,n]-Subscript[t, 1,n]*E^(I*Subscript[k1, z,n]*tc))*Subscript[b1, n,m],{n,-j,j}]-Sum[Subscript[k2, z,n]*(-Subscript[t, 2,n]+Subscript[r, 2,n]*E^(I*Subscript[k2, z,n]*(t-tc)))*Subscript[b2, n,m],{n,-j,j}]==0,{m,-j,j}],

Table[Sum[(Subscript[r, 2,n]+Subscript[t, 2,n]*E^(I*Subscript[k2, z,n]*(t-tc)))*Subscript[b2, n,m],{n,-j,j}]-Sum[(Subscript[t, 3,n]+Subscript[r, 3,n]*E^(I*Subscript[k3, z,n]*tc))*Subscript[b3, n,m],{n,-j,j}]==0,{m,-j,j}],Table[Sum[Subscript[k2, z,n]*(Subscript[r, 2,n]-Subscript[t, 2,n]*E^(I*Subscript[k2, z,n]*(t-tc)))*Subscript[b2, n,m],{n,-j,j}]-Sum[Subscript[k3, z,n]*(-Subscript[t, 3,n]+Subscript[r, 3,n]*E^(I*Subscript[k3, z,n]*tc))*Subscript[b3, n,m],{n,-j,j}]==0,{m,-j,j}],


Table[Sum[(Subscript[t, 3,n]*E^(I*Subscript[k3, z,n]*tc)+Subscript[r, 3,n])*Subscript[b3, n,m],{n,-j,j}]-Subscript[t, 4,m]==0,{m,-j,j}],Table[Sum[Subscript[k3, z,n]*(-Subscript[t, 3,n]*E^(I*Subscript[k3, z,n]*tc)+Subscript[r, 3,n])*Subscript[b3, n,m],{n,-j,j}]+Subscript[p, m]*Subscript[t, 4,m]==0,{m,-j,j}]};

arr=Flatten[Arr];

lis={Table[Subscript[r, 0,n],{n,-j,j}],Table[Subscript[r, 1,n],{n,-j,j}],Table[Subscript[t, 1,n],{n,-j,j}],Table[Subscript[r, 2,n],{n,-j,j}],Table[Subscript[t, 2,n],{n,-j,j}],Table[Subscript[r, 3,n],{n,-j,j}],Table[Subscript[t, 3,n],{n,-j,j}],Table[Subscript[t, 4,n],{n,-j,j}] };
list=Flatten[lis];

Sol=Solve[arr,list];
Table[Subscript[RF, s,h,n]=Subscript[R, 0,n]=Sol[[1]][[n+j+1,2]],{n,-j,j}];
Table[Subscript[R, 1,n]=Sol[[1]][[n+j+1+(2j+1),2]],{n,-j,j}];
Table[Subscript[T, 1,n]=Sol[[1]][[n+j+1+2*(2j+1),2]],{n,-j,j}];
Table[Subscript[R, 2,n]=Sol[[1]][[n+j+1+3*(2j+1),2]],{n,-j,j}];
Table[Subscript[T, 2,n]=Sol[[1]][[n+j+1+4*(2j+1),2]],{n,-j,j}];
Table[Subscript[R, 3,n]=Sol[[1]][[n+j+1+5*(2j+1),2]],{n,-j,j}];
Table[Subscript[T, 3,n]=Sol[[1]][[n+j+1+6*(2j+1),2]],{n,-j,j}];
Table[Subscript[TF, s,h,n]=Subscript[T, 2,n]=Sol[[1]][[n+j+1+7*(2j+1),2]],{n,-j,j}];
(*Condition to avoid unphysical and decaying solutions:*)
Table[If[Im[Subscript[X, s,h,n]]!=0,Subscript[Rf, s,h,n]=0,Subscript[Rf, s,h,n]=Subscript[RF, s,h,n]],{n,-j,j}];
Table[If[Im[Subscript[Pt, s,h,n]]!=0,Subscript[Tf, s,h,n]=0,Subscript[Tf, s,h,n]=Subscript[TF, s,h,n]],{n,-j,j}];
(*Print[Subscript[k, y,s,1],Subscript[k, z,s,1]];*)
(*Print[s];*)
s=s+1];
Print[h];
h=h+1]


For[i=-15, i<0, i++,bigmatrxw=Table[{Re[Subscript[X, s,1,i]],((Re[Subscript[X,s,1,i]]/Re[Subscript[k, z,s,1]])*Abs[Subscript[Rf, s,1,i]])^2},{s,1,443}];
newdata=MovingAverage[bigmatrxw,3];
Export["/Users/hadirahmaninejad/Desktop/x-ray-main/Project of x-rays/Fits/Optimization/multiple_condition_simultaneous/dyt/rod="<>ToString[-i]<>".txt",newdata,"Table"];Print[i]]



(*R=Import["/Users/hadirahmaninejad/Desktop/Ashkar's lab projects/Optimization/Process_of_fitting/Main/test/three_variables/samplenewly.txt","Table"][[All,{1,2}]];
Show[ListPlot[R,PlotStyle\[Rule] {Gray,Solid},PlotLegends\[Rule]{"\[Rho]f>\[Rho]s"}],
AxesLabel\[Rule]{Style[Wave number [1/nm],FontSize\[Rule] 14, FontColor\[Rule]Black],Style["R"^2,FontSize\[Rule] 14, FontColor\[Rule]Black]}, AxesStyle \[Rule] 10]*)


(*R=Import["/Users/hadirahmaninejad/Desktop/Ashkar's lab projects/Optimization/Process_of_fitting/Main/test/three_variables/sample.txt","Table"][[All,{1,2}]];

Show[ListPlot[R,ScalingFunctions\[Rule]"Log",PlotStyle\[Rule] {Gray,Solid},PlotRange\[Rule] {{0,0.01},{0,0.1}},PlotLegends\[Rule]{"\[Rho]f>\[Rho]s"}],PlotRange\[Rule]{{0,0.03},All},AxesLabel\[Rule]{Style[Wave number [1/nm],FontSize\[Rule] 14, FontColor\[Rule]Black],Style["R"^2,FontSize\[Rule] 14, FontColor\[Rule]Black]}, AxesStyle \[Rule] 10,AxesOrigin \[Rule] {0,-18}]*)

