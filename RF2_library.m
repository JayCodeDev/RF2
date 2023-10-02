%RF2 library of funcitons

function ExamPrep()

fprintf("\n***************\nQuestion #2\n***************\n");
Qn=4;
[C1,C2,R]=Qfactor(Qn);
figure(1)
SmithCircle(real(C1),imag(C1),R);
SmithCircle(real(C2),imag(C2),R);

fprintf("\n***************\nQuestion #3\n***************\n");
S11=Pol2Rectd(0.44,-74);
S12=Pol2Rectd(0.04,64); %Unilateral when S12=0
S21=Pol2Rectd(23,126);
S22=Pol2Rectd(0.67,-45);
GammaL=Pol2Rectd(.5,0);
GammaIN=S11+((S12)*(S21)*GammaL)/(1-(S22)*GammaL);
GammaS=conj(S11+((S12)*(S21)*GammaL)/(1-(S22)*GammaL))
[A,Theta] = Rect2Pold(real(GammaS),imag(GammaS))

[Del,K]=Del_K(S11,S12,S21,S22);
[B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);

fprintf("\n***************\nQuestion #4\n***************\n");
fprintf("\n***************\nPart a)\n***************\n");
S11 = Pol2Rectd(0.5,-70);
S12 = Pol2Rectd(0.2,-10);
S21 = Pol2Rectd(5,80);
S22 = Pol2Rectd(0.1,-30);
GPorAdB=12;
Wanted_Circle="In_A"; %"Out_P" or "In_A"
[Del, K]=Del_K(S11,S12,S21,S22);
fprintf("\nGAdB = %f",GPorAdB);
[C,R] = BiLat_GainCircles(S11,S12,S21,S22,GPorAdB,Wanted_Circle,Del,K);
figure(2)
SmithCircle(real(C),imag(C),R);

GammaOPT = Pol2Rectd(0.3,-135);
FmindB = 1.2;
RnDenorm = 25;
FidB = 1.5;
fprintf("\nFidB = %f",FidB);
Z0 = 50;
[C,R] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(C),imag(C),R);
FidB = 2;
fprintf("\nFidB = %f",FidB);
[C,R] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(C),imag(C),R);

% GammaS = -0.3-1i*.41
GammaS = Pol2Rectd(0.50961,-105.945);
fprintf("\n***************\nPart b)\n***************\n");
S11or22 = S22;
S22or11 = S11;
GammaLorS = GammaS;
InOut="OUT"; %"IN" or "OUT"
[GammaOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);
GammaINorOut=GammaOUT;
GammaSorL=conj(GammaINorOut);
[GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL);
[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB);

S11or22 = S11;
S22or11 = S22;
GammaLorS = GammaL;
InOut="IN"; %"IN" or "OUT"
[GammaIN] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);
GammaINorOut=GammaIN;
GammaSorL=GammaS;
[GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL);
[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB);
SmithCircle(real(VSWRINorOUTcent),imag(VSWRINorOUTcent),VSWRINorOUTrad);
GammaS
GammaL
GammaIN
GammaOUT
VSWRIN=VSWRINorOUT

fprintf("\n***************\nQuestion #5\n***************\n");

end
function [C1,C2,R]=Qfactor(Qn)
% Qn=4;
% [C1,C2,R]=Qfactor(Qn);

C=1/Qn;
C1=1i*C;
C2=-1i*C;
R=sqrt(1+1/(Qn^2));

fprintf("\n*** Qn Setup ***\nWith Qn = %f\nCenter is +/- %f\nRadius is %f\n",Qn,C,R);

end
function [Gamma] = Z2Gamma(ZL,Z0)

Gamma=(ZL-Z0)/(ZL+Z0);

end
function [ZL] = Gamma2Z(Gamma,Z0)

ZL=Z0*(1+Gamma)/(1-Gamma);

end
function [F,M1,M2]=Fcalc2Amps(G1Lin,F1Lin,G2Lin,F2Lin)
% G1Lin = 3;  %GTU Max Value
% F1Lin = 2;
% G2Lin = 3;  %GTU Max Value
% F2Lin = 2.5;
% [F,M1,M2]=Fcalc2Amps(G1Lin,F1Lin,G2Lin,F2Lin)

F = F1Lin+(F2Lin-1)/(G1Lin);
M1 = (F1Lin-1)/(1-1/G1Lin);
M2 = (F2Lin-1)/(1-1/G2Lin);

fprintf("\nThe Amplifiers have the following M values\nKeep in mind the lowest should go first\n");
fprintf("Amp 1: M1 = %f\n",M1);
fprintf("Amp 2: M2 = %f\n",M2);
fprintf("\nThe current setup has a noise figure F = %f\n",F);

end
function [F,M1,M2,M3]=Fcalc3Amps(G1Lin,F1Lin,G2Lin,F2Lin,G3Lin,F3Lin)
% G1Lin = 3;  %GTU Max Value
% F1Lin = 2;
% G2Lin = 3;  %GTU Max Value
% F2Lin = 2.5;
% G3Lin = 3;  %GTU Max Value
% F3Lin = 1;
% [F,M1,M2,M3]=Fcalc3Amps(G1Lin,F1Lin,G2Lin,F2Lin,G3Lin,F3Lin)

F=F1Lin+(F2Lin-1)/(G1Lin)+(F3Lin-1)/(G1Lin*G2Lin);
M1 = (F1Lin-1)/(1-1/G1Lin);
M2 = (F2Lin-1)/(1-1/G2Lin);
M3 = (F3Lin-1)/(1-1/G3Lin);

fprintf("\nThe Amplifiers have the following M values\nKeep in mind the lowest should go first\n");
fprintf("Amp 1: M1 = %f\n",M1);
fprintf("Amp 2: M2 = %f\n",M2);
fprintf("Amp 3: M2 = %f\n",M3);
fprintf("\nThe current setup has a noise figure F = %f\n",F);

end
function [F,M1,M2,M3,M4]=Fcalc4Amps(G1Lin,F1Lin,G2Lin,F2Lin,G3Lin,F3Lin,G4Lin,F4Lin)
% G1Lin = 3;  %GTU Max Value
% F1Lin = 2;
% G2Lin = 3;  %GTU Max Value
% F2Lin = 2.5;
% G3Lin = 3;  %GTU Max Value
% F3Lin = 1;
% G4Lin = 3;  %GTU Max Value
% F4Lin = 3;
% [F]=Fcalc4Amps(G1Lin,F1Lin,G2Lin,F2Lin,G3Lin,F3Lin,G4Lin,F4Lin)

F=F1Lin+(F2Lin-1)/(G1Lin)+(F3Lin-1)/(G1Lin*G2Lin)+(F4Lin-1)/(G1Lin*G2Lin*G3Lin);
M1 = (F1Lin-1)/(1-1/G1Lin);
M2 = (F2Lin-1)/(1-1/G2Lin);
M3 = (F3Lin-1)/(1-1/G3Lin);
M4 = (F4Lin-1)/(1-1/G4Lin);

fprintf("\nThe Amplifiers have the following M values\nKeep in mind the lowest should go first\n");
fprintf("Amp 1: M1 = %f\n",M1);
fprintf("Amp 2: M2 = %f\n",M2);
fprintf("Amp 3: M2 = %f\n",M3);
fprintf("Amp 4: M2 = %f\n",M4);
fprintf("\nThe current setup has a noise figure F = %f\n",F);

end
function Quiz07()

fprintf("\n***************\nPart A)\n***************\n");
S11=Pol2Rectd(0.6,-170);
S12=Pol2Rectd(0.05,50); 
S21=Pol2Rectd(8,80);
S22=Pol2Rectd(0.45,-10);
InOut="Out"; %"In" or "Out"
[C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut)
figure(1)
SmithCircle(real(C),imag(C),R);

InOut="In"; %"In" or "Out"
[C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut)
figure(2)
SmithCircle(real(C),imag(C),R);

fprintf("\n***************\nPart B)\n***************\n");
[G_MSG] = MaxStabGain(S12,S21)

fprintf("\n***************\nPart C)\n***************\n");
GPorAdB=19;
Wanted_Circle="Out_P"; %"Out_P" or "In_A"
[Del, K]=Del_K(S11,S12,S21,S22);
[C,R] = BiLat_GainCircles(S11,S12,S21,S22,GPorAdB,Wanted_Circle,Del,K);
SmithCircle(real(C),imag(C),R);

fprintf("\n***************\nPart D)\n***************\n");
[AMatch,ThetaMatch,Gammas,GammaSorLoptions] = MatchingCircle(real(C),imag(C),R);
GammaL=Pol2Rectd(AMatch,ThetaMatch)

fprintf("\n***************\nPart E)\n***************\n");
S11or22 = S11;
S22or11 = S22;
GammaLorS = GammaL;
InOut="IN"; %"IN" or "OUT"
[GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);

fprintf("\n***************\nPart F)\n***************\n");
GammaS=conj(GammaINorOUT)+0.2

fprintf("\n***************\nPart G)\n***************\n");
[A,Theta] = Rect2Pold(real(GammaS),imag(GammaS))

fprintf("\n***************\nPart H)\n***************\n");
S11or22 = S22;
S22or11 = S11;
GammaLorS = GammaS;
InOut="OUT"; %"IN" or "OUT"
[GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);

fprintf("\n***************\nPart I)\n***************\n");
GammaINorOut=GammaINorOUT;
GammaSorL=GammaL;
[GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL)
[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB)
figure(1)
C=VSWRINorOUTcent;
R=VSWRINorOUTrad;
SmithCircle(real(C),imag(C),R);

fprintf("\n***************\nPart J)\n***************\n");
fprintf("\nWe can change GammaS to give us a large VSWR_IN\n");

end
function [C,R] = BiLat_GainCircles(S11,S12,S21,S22,GPorAdB,Wanted_Circle,Del,K)
% S11 = Pol2Rectd(0.6,-140);
% S12 = Pol2Rectd(0,0);
% S21 = Pol2Rectd(4,-90);
% S22 = Pol2Rectd(0.4,-60);
% GPorAdB=17;
% Wanted_Circle="Out_P"; %"Out_P" or "In_A"
% [Del, K]=Del_K(S11,S12,S21,S22);
% [C,R] = BiLat_GainCircles(S11,S12,S21,S22,GPorAdB,Wanted_Circle,Del,K);
% SmithCircle(real(C),imag(C),R);

B1 = 1+norm(S11)^2-norm(S22)^2-norm(Del)^2;
B2 = 1+norm(S22)^2-norm(S11)^2-norm(Del)^2;
C1 = S11-Del*conj(S22);
C2 = S22-Del*conj(S11);
GLin = 10^(GPorAdB/10);
g=GLin/(norm(S21)^2);
if Wanted_Circle == "Out_P"
 C = (g*conj(C2))/(1+g*(norm(S22)^2-norm(Del)^2));
 R = sqrt(1-2*K*norm(S12*S21)*g+norm(S12*S21)^2*g^2)/(norm(1+g*(norm(S22)^2-norm(Del)^2)));
elseif Wanted_Circle == "In_A"
 C = (g*conj(C1))/(1+g*(norm(S11)^2-norm(Del)^2));
 R = sqrt(1-2*K*norm(S12*S21)*g+norm(S12*S21)^2*g^2)/norm(1+g*(norm(S11)^2-norm(Del)^2));
end
[A,Theta] = Rect2Pold(real(C),imag(C));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(C),imag(C),R)
end
function [G_MSG] = MaxStabGain(S12,S21)
% S12 = Pol2Rectd(0,0);
% S21 = Pol2Rectd(4,-90);
% [G_MSG] = MaxStabGain(S12,S21)

G_MSG=norm(S21)/norm(S12);

end
function RF2HW9()

fprintf("\n*** Problem 1) ***\n");

%From Example 4.3.4
S11 = Pol2Rectd(0.7,-105);
S12 = Pol2Rectd(0.11,20);
S21 = Pol2Rectd(3,75);
S22 = Pol2Rectd(0.46,-70);
FmindB = 0.8; %dB
GammaOPT = Pol2Rectd(0.7,55);
RnDenorm = 0.95*50; %Ohm
VDS = 3;    %V
IDS = 0.020;    %A
f=6e9;  %GHz

%Assume
Z0=50;

%Target
VSWRIN=2;

[Del K]=Del_K(S11,S12,S21,S22);

InOut="In"; %"In" or "Out"
[C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut);
SmithCircle(real(C),imag(C),R);

% InOut="Out"; %"In" or "Out"
% [C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut);
% SmithCircle(real(C),imag(C),R);

GMSGLin=norm(S21)/norm(S12);
GMSGdB=10*log10(GMSGLin);
fprintf("\nGMSG = %f (or %f dB)\n",GMSGLin,GMSGdB);

C1=S11-Del*conj(S22);
C2=S22-Del*conj(S11);
Ci = C2;
Sii = S22;

fprintf("\nGP=10 dB Circle");
GPorAdB=10;
[C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);
SmithCircle(real(C),imag(C),R);

fprintf("\nGA=10 dB Circle Converted");
Type='out';
[C, R] = Convert_Circles(S11, S12, S21, S22, C, R, Type);
SmithCircle(real(C),imag(C),R);
[A,Theta] = Rect2Pold(real(C),imag(C));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(C),imag(C),R)

fprintf("\nGP=12 dB Circle");
GPorAdB=12;
[C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);
SmithCircle(real(C),imag(C),R);

fprintf("\nGA=12 dB Circle Converted");
Type='out';
[C, R] = Convert_Circles(S11, S12, S21, S22, C, R, Type);
SmithCircle(real(C),imag(C),R);
[A,Theta] = Rect2Pold(real(C),imag(C));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(C),imag(C),R)

fprintf("\nF=1.2 dB Circle");
FidB = 1.2;
[C,R] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(C),imag(C),R);

fprintf("\nF=1.5 dB Circle");
FidB = 1.5;
[C,R] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(C),imag(C),R);

fprintf("\nVSWRIN=2 Circle");
VSWRINorOUT = VSWRIN;
[GammaAorB] = VSWRINorOUTtoGammaAorB(VSWRINorOUT);

S11or22 = S11;
S22or11 = S22;
GammaLorS = Pol2Rectd(0.42,-57.3);
InOut="IN"; %"IN" or "OUT"
[GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);

[VSWRINorOUT,C,R] = VSWRINorOUT_CR_Calc(GammaINorOUT,GammaAorB);
SmithCircle(real(C),imag(C),R);

fprintf("\n*** Problem 2) ***\n");

F1 =10^(0.5 /10); %Linear
GA1=10^(14  /10); %Linear
F2 =10^(0.9 /10); %Linear
GA2=10^(14  /10); %Linear
F3 =10^(1.1 /10); %Linear
GA3=10^(14  /10); %Linear
F4 =10^(7.5 /10); %Linear
GA4=10^(-7.5/10);   %Linear
F5 =10^(4.5 /10); %Linear
GA5=10^(16  /10); %Linear
F6 =10^(4.5 /10); %Linear
GA6=10^(16  /10); %Linear

FLin=F1 + (F2-1)/(GA1) + (F3-1)/(GA1*GA2) + (F4-1)/(GA1*GA2*GA3) + (F5-1)/(GA1*GA2*GA3*GA4) + (F6-1)/(GA1*GA2*GA3*GA4*GA5);
FdB=10*log10(FLin);
GLin=GA1*GA2*GA3*GA4*GA5*GA6;
GdB=10*log10(GLin);

fprintf("\nF = %f (or %f dB)\n",FLin,FdB);
fprintf("G = %f (or %f dB)\n",GLin,GdB);

fprintf("\n*** Problem 3) ***\n");

S11 = Pol2Rectd(0.97,-8);
S21 = Pol2Rectd(7.7,177);
S22 = Pol2Rectd(0.97,-7);

GT=10; %dB
Z0=50;

S21desired=10^(GT/20);
R2=Z0*(1+S21desired);
fprintf("\nR2 = %f\n",R2);

fprintf("\n*** Problem 4) ***\n");
fprintf("a)\n");

C = 100e-12;    %F
R = 5;  %Ohms
fa = 400e6;
fb = 600e6;
f0 = (fa+fb)/2;

Q1=1/(2*pi*f0*R*C);
Q2=2*pi*f0/(2*pi*fb-2*pi*fa);
GammaX=exp(-pi*Q2/Q1);
fprintf("\nQ1 = %f\nQ2 = %f\nGammaX = %fe-6\n",Q1,Q2,GammaX*1e6);

fprintf("\nb)\n");

C = 1e-12;    %F
R = 50;  %Ohms
fa = 6e9;
fb = 12e9;
f0 = (fa+fb)/2;

XC=1/(2*pi*f0*C);
Q1=R/XC;
Q2=2*pi*f0/(2*pi*fb-2*pi*fa);
GammaX=exp(-pi*Q2/Q1);
fprintf("\nQ1 = %f\nQ2 = %f\nGammaX = %f\n",Q1,Q2,GammaX);

end
function [Center, Radius] = Convert_Circles(S11, S12, S21, S22, C, r, Type)
Del = (S11)*(S22) - (S12)*(S21);
MagDel = norm(Del);
Kay = (1-(norm(S11)^2)-(norm(S22)^2)+MagDel^2)/(2*norm((S12)*(S21)));

if Type == 'out'
    Centeri = ((1-S22*C)*conj(S11-Del*C)-r^2*conj(Del)*S22)/abs((abs(1-S22*C)^2-r^2*abs(S22)^2));
    ri = (r*abs(S12*S21))/abs(abs(1-S22*C)^2-r^2*abs(S22)^2);
    Center = Centeri;
    Radius = ri;
elseif Type == 'in'
    Centero = ((1-S11*C)*conj(S22-Del*C)-r^2*conj(Del)*S11)/abs((abs(1-S11*C)^2-r^2*abs(S11)^2));
    ro = (r*abs(S12*S21))/abs(abs(1-S11*C)^2-r^2*abs(S11)^2);
    Radius = ro;
end
end
function RF2HW8()

S11 = Pol2Rectd(0.6,-170);
S12 = Pol2Rectd(0.05,16);
S21 = Pol2Rectd(2,30);
S22 = Pol2Rectd(0.5,-95);
FmindB = 2.5; %dB
GammaOPT = Pol2Rectd(0.5,145);
RnDenorm = 5; %Ohm
Z0=50;

[Del, K]=Del_K(S11,S12,S21,S22);
[B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);

GPorAdB = GTmax-3;
Ci = C2;
Sii = S22;
[C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);

SmithCircle(real(C),imag(C),R);

fprintf("\nFi = 3 dB\n");

FidB=3;
[CFi,RFi] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(CFi),imag(CFi),RFi);

fprintf("\nFi = 4 dB\n");

FidB=4;
[CFi,RFi] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);
SmithCircle(real(CFi),imag(CFi),RFi);

[FdB,FLin] = NoiseFigureCalc(GammaMS,GammaOPT,FmindB,RnDenorm,Z0);

end
function [FdB,FLin] = NoiseFigureCalc(GammaMS,GammaOPT,FmindB,RnDenorm,Z0)
% GammaMS = Pol2Rectd(0.5,145);
% GammaOPT = Pol2Rectd(0.5,145);
% FmindB = 2.5;
% RnDenorm = 5;
% FidB = 3;
% Z0 = 50;
% [F] = NoiseFigureCalc(GammaMS,GammaOPT,FmindB,RnDenorm,Z0);

FminLin = 10^(FmindB/10);
RnNorm = RnDenorm/Z0;

gs = (1-GammaMS)/(1+GammaMS);
gopt = (1-GammaOPT)/(1+GammaOPT);

F = FminLin+RnNorm/gs*abs(gs-gopt)^2;
FLin = F;
FdB = 10*log10(F);

fprintf("\nF = %f + j%f\n|F| = %f Linear (%f dB)\n",real(F),imag(F),norm(F),FdB);

end
function [CFi,RFi] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0)
% GammaOPT = Pol2Rectd(0.5,145);
% FmindB = 2.5;
% RnDenorm = 5;
% FidB = 3;
% Z0 = 50;
% [CFi,RFi] = NoiseCR(GammaOPT,FmindB,RnDenorm,FidB,Z0);

FiLin=10^(FidB/10);
FminLin=10^(FmindB/10);
RnNorm=RnDenorm/Z0;

Ni=((FiLin-FminLin)*norm(1+GammaOPT)^2)/(4*RnNorm);
CFi=GammaOPT/(1+Ni);
RFi=sqrt((Ni^2)+Ni*(1-norm(GammaOPT)^2))/(1+Ni);

[A,Theta] = Rect2Pold(real(CFi),imag(CFi));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(CFi),imag(CFi),RFi)
end
function RF2HW7()
syms x y 
figcount=0;
% 
fprintf("\n******************\n*** Problem 1) ***\n******************\n");
S11=Pol2Rectd(0.7,-65);
S12=Pol2Rectd(0.03,60); 
S21=Pol2Rectd(3.2,110);
S22=Pol2Rectd(0.8,-30);
f=2e9;
w=2*pi*f;
Gp=10;

[Del,K]=Del_K(S11,S12,S21,S22);
[B1,B2,C1,C2,GTmax,GTmin,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);

Ci = C2;
Sii = S22;
GPorAdB=Gp;
[C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);

figcount=figcount+1;
figure(figcount)
title('Gain Circle');
[A,Theta,Gammas] = MatchingCircle(real(C),imag(C),R);

GammaL=Pol2Rectd(A,Theta);
[GammaIN] = UniGammaINCalc(S11,S12,S21,S22,GammaL);
GammaS=conj(GammaIN);
[A,Theta] = Rect2Pold(real(GammaL),imag(GammaL));

GammaINorOut=GammaIN;
GammaSorL=GammaS;
[GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL);

[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB);
VSWRIN=VSWRINorOUT;

figcount=figcount+1;
figure(figcount)
title('VSWRIN Circle');
SmithCircle(real(VSWRINorOUTcent),imag(VSWRINorOUTcent),VSWRINorOUTrad);

GammaOUT=S22+(S12*S21*GammaS)/(1-S11*GammaS);
GammaINorOut=GammaOUT;
GammaSorL=GammaL;
[GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL);

[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB);
VSWROUT=VSWRINorOUT;

figcount=figcount+1;
figure(figcount)
title('VSWROUT Circle');
SmithCircle(real(VSWRINorOUTcent),imag(VSWRINorOUTcent),VSWRINorOUTrad);

fprintf("\nVSWRIN = %f\nVSWROUT = %f\n",VSWRIN,VSWROUT);

fprintf("\n******************\n*** Problem 2) ***\n******************\n");
f=12e9;     %Hz
w=2*pi*f;   %Rad/S
VDS=3.5;    %V
ID=25e-3;   %A
S11=Pol2Rectd(0.6,36);
S12=Pol2Rectd(0.14,-85); 
S21=Pol2Rectd(2.3,-80);
S22=Pol2Rectd(0.15,45);

fprintf("\n******************\n*** Part a)    ***\n******************\n");

[Del,K]=Del_K(S11,S12,S21,S22);
[B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);

Gp=GTmaxdB-1

Ci = C2;
Sii = S22;
GPorAdB=Gp;
[C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);

fprintf("\n******************\n*** Part b)    ***\n******************\n");

figcount=figcount+1;
figure(figcount)
title('Gain Circle');
[AMatch,ThetaMatch,Gammas,GammaSorLoptions] = MatchingCircle(real(C),imag(C),R);
GammaL=GammaSorLoptions

S11or22 = S11;
S22or11 = S22;
GammaLorS = GammaL;
InOut="IN"; %"IN" or "OUT"
[GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);
GammaIN=GammaINorOUT

VSWRINorOUT=1.8;
[GammaAorB] = VSWRINorOUTtoGammaAorB(VSWRINorOUT);
GammaA=GammaAorB

GammaINorOut=GammaIN;
[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaIN,GammaA);

figcount=figcount+1;
figure(figcount)
title('VSWR Circles');
SmithCircle(real(VSWRINorOUTcent(1,1)),imag(VSWRINorOUTcent(1,1)),VSWRINorOUTrad(1,1));
SmithCircle(real(VSWRINorOUTcent(2,1)),imag(VSWRINorOUTcent(2,1)),VSWRINorOUTrad(2,1));
SmithCircle(real(VSWRINorOUTcent(3,1)),imag(VSWRINorOUTcent(3,1)),VSWRINorOUTrad(3,1));
SmithCircle(real(VSWRINorOUTcent(4,1)),imag(VSWRINorOUTcent(4,1)),VSWRINorOUTrad(4,1));

[xout1,yout1] = circcirc(real(C),imag(C),R,real(VSWRINorOUTcent(1,1)),imag(VSWRINorOUTcent(1,1)),VSWRINorOUTrad(1,1));
[xout2,yout2] = circcirc(real(C),imag(C),R,real(VSWRINorOUTcent(2,1)),imag(VSWRINorOUTcent(2,1)),VSWRINorOUTrad(2,1));
[xout3,yout3] = circcirc(real(C),imag(C),R,real(VSWRINorOUTcent(3,1)),imag(VSWRINorOUTcent(3,1)),VSWRINorOUTrad(3,1));
[xout4,yout4] = circcirc(real(C),imag(C),R,real(VSWRINorOUTcent(4,1)),imag(VSWRINorOUTcent(4,1)),VSWRINorOUTrad(4,1));

GammaSOptions = [xout1' yout1';xout2' yout2';xout3' yout3';xout4' yout4'];

[A1,Theta1] = Rect2Pold(GammaSOptions(1,1),GammaSOptions(1,2));
fprintf("GammaS Option 1: %f @ %f\n",A1,Theta1);
[A2,Theta2] = Rect2Pold(GammaSOptions(2,1),GammaSOptions(2,2));
fprintf("GammaS Option 2: %f @ %f\n",A2,Theta2);

fprintf("\n******************\n*** Part C)    ***\n******************\n");

S11or22 = S22;
S22or11 = S11;
GammaLorS = [Pol2Rectd(A1,Theta1);Pol2Rectd(A2,Theta2)];
InOut="OUT"; %"IN" or "OUT"
[GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);
GammaOUT=GammaINorOUT

[GammaAorB] = GammaAorBCalc(GammaOUT,GammaL(1,1));
GammaB=GammaAorB

[VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaOUT,GammaB);
VSWROUT=VSWRINorOUT

fprintf("\n******************\n*** Problem 3) ***\n******************\n");
VCC = 12;   %V
IC = 20e-3; %A
VCE = 5;    %V
Beta = 125;
VBE = 0.75; %V

VC = 9 %V    %Arbitrary pick

IB = IC / Beta
I1 = IB + IC
R1 = (VCC - VC)/I1
R2 = (VCE - VBE)/IB
R3 = (VCC - I1*R1-VCE)/I1


end
function [GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut)
% S11or22 = Pol2Rectd(0.7,-65);
% S22or11 = Pol2Rectd(0.8,-30);
% S12 = Pol2Rectd(0.03,60);
% S21 = Pol2Rectd(3.2,110);
% GammaLorS = Pol2Rectd(0.5,40);
% InOut="IN"; %"IN" or "OUT"
% [GammaINorOUT] = GammaInorOut(S11or22,S22or11,S12,S21,GammaLorS,InOut);

GammaINorOUT=S11or22+(((S12*S21).*GammaLorS)./(1-S22or11.*GammaLorS));

ArrayCheck=length(GammaLorS)>1;
fprintf("\n");
if ArrayCheck == 0
    [A,Theta] = Rect2Pold(real(GammaINorOUT),imag(GammaINorOUT));
    fprintf("Gamma%s = %f @ %f\n",InOut,A,Theta);
elseif ArrayCheck == 1
    count=1;
    while count <= length(GammaLorS)
        
        [A,Theta] = Rect2Pold(real(GammaINorOUT(count,1)),imag(GammaINorOUT(count,1)));
        fprintf("Gamma%s option %i = %f @ %f\n",InOut,count,A,Theta);
        
        count=count+1;
    end
end

end
function [GammaAorB] = VSWRINorOUTtoGammaAorB(VSWRINorOUT)

GammaAorB=(VSWRINorOUT-1)/(VSWRINorOUT+1);

end
function [C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB)
% S11=Pol2Rectd(0.7,-65);
% S12=Pol2Rectd(0.03,60); 
% S21=Pol2Rectd(3.2,110);
% S22=Pol2Rectd(0.8,-30);
% GPorAdB = 10;
% [Del,K]=Del_K(S11,S12,S21,S22);
% [B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);
% Ci = C2;
% Sii = S22;
% GPorAdB=Gp;
% [C,R,gPorA,GPorAlin] = Bilateral_CR_Calc(Sii,S12,S21,Ci,K,Del,GPorAdB);

GPorAlin=10^(GPorAdB/10);
gPorA=GPorAlin/norm(S21)^2;

C=gPorA*conj(Ci)/(1+gPorA*(norm(Sii)^2-norm(Del)^2));
R=sqrt(1-2*K*norm(S12*S21)*gPorA+norm(S12*S21)^2*gPorA^2)/norm(1+gPorA*(norm(Sii)^2-norm(Del)^2));

[A,Theta] = Rect2Pold(real(C),imag(C));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(C),imag(C),R)
end
function [B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del)
% S11=Pol2Rectd(0.7,-65);
% S12=Pol2Rectd(0.03,60); 
% S21=Pol2Rectd(3.2,110);
% S22=Pol2Rectd(0.8,-30);
% [Del,K]=Del_K(S11,S12,S21,S22);
% [B1,B2,C1,C2,GTmax,GTmaxdB,GTmin,GTmindB,GammaMS,GammaML] = GammaMSorMLCalc(S11,S12,S21,S22,K,Del);

posneg=1;

B1=1+norm(S11)^2-norm(S22)^2-norm(Del)^2;
B2=1+norm(S22)^2-norm(S11)^2-norm(Del)^2;
C1=S11-Del*conj(S22);
C2=S22-Del*conj(S11);
fprintf("\n");
if K>1 && norm(Del)<1
    posneg=-1;
    GTmax = norm(S21)/norm(S12)*(K-sqrt((K^2)-1));
    GTmaxdB=10*log10(GTmax);
    GTmin = 0;
    GTmindB=10*log10(GTmin);
    fprintf("GT,max = GP,max = GA,max = %f (or %f dB)\nGT,min = %f (or %f dB)\n",GTmax,GTmaxdB,GTmin,GTmindB);
elseif K>1 && norm(Del)>1
    posneg=1;
    GTmax = inf;
    GTmaxdB=10*log10(GTmax);
    GTmin = norm(S21)/norm(S12)*(K+sqrt((K^2)-1));
    GTmindB=10*log10(GTmin);
    fprintf("GT,max = %f (or %f dB)\nGT,min = %f (or %f dB)",GTmax,GTmaxdB,GTmin,GTmindB);
elseif K==1
    GTmax=norm(S21)/norm(S12);
    GTmaxdB=10*log10(GTmax);
    GTmin = 0;
    GTmindB=10*log10(GTmin);
    fprintf("Max stable gain with K=1 is GMSG = %f (or %f dB)\nGT,min = %f (or %f dB)\n",GTmax,GTmaxdB,GTmin,GTmindB);
elseif K<1
    GTmax="NA";
    GTmaxdB="NA";
    GTmin ="NA";
    GTmindB="NA";
    fprintf("Simultaneous match is NOT possible with |GammaML|<1 and |GammaMS|<1\n");
end
GammaMS=(B1+posneg*sqrt(B1^2-4*norm(C1)^2))/(2*C1);
GammaML=(B2+posneg*sqrt(B2^2-4*norm(C2)^2))/(2*C2);

[A,Theta] = Rect2Pold(real(GammaMS),imag(GammaMS));
fprintf("GammaMS = %f @ %f\n",A,Theta);
[A,Theta] = Rect2Pold(real(GammaML),imag(GammaML));
fprintf("GammaML = %f @ %f\n",A,Theta);
end
function [VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB)
% GammaINorOut=GammaIN;
% GammaSorL=GammaS;
% [GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL)
% [VSWRINorOUT,VSWRINorOUTcent,VSWRINorOUTrad] = VSWRINorOUT_CR_Calc(GammaINorOut,GammaAorB)

VSWRINorOUT = (1+abs(GammaAorB))./(1-abs(GammaAorB));
VSWRINorOUTcent = (conj(GammaINorOut).*(1-abs(GammaAorB).^2))./(1-abs(GammaAorB.*GammaINorOut).^2);
VSWRINorOUTrad = (abs(GammaAorB).*(1-abs(GammaINorOut).^2))./(1-abs(GammaAorB.*GammaINorOut).^2);

[A,Theta] = Rect2Pold(real(VSWRINorOUTcent),imag(VSWRINorOUTcent));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(VSWRINorOUTcent),imag(VSWRINorOUTcent),VSWRINorOUTrad)
end
function [GammaAorB] = GammaAorBCalc(GammaINorOut,GammaSorL)
% GammaINorOut=GammaIN;
% GammaSorL=GammaS;
% [GammaAorBCalc] = GammaAorBCalc(GammaINorOut,GammaSorL)

GammaAorB=abs((GammaINorOut-conj(GammaSorL))./(1-GammaINorOut.*GammaSorL));

end
function [h] = SmithCircle(x,y,r)
% SmithCircle(real(C),imag(C),R);
hold on
XMIN=-1;
XMAX=1;
YMIN=-1;
YMAX=1;
axis([XMIN XMAX YMIN YMAX]);
th = 0:pi/50:2*pi;
xunit1 = 1 * cos(th) + 0;
yunit1 = 1 * sin(th) + 0;
h = plot(xunit1, yunit1,'r');
xunit2 = 0.5 * cos(th) + 0.5;
yunit2 = 0.5 * sin(th) + 0;
h = plot(xunit2, yunit2,'r');
xunit3 = 0.5 * cos(th) - 0.5;
yunit3 = 0.5 * sin(th) + 0;
h = plot(xunit3, yunit3,'b');
xunit4 = r * cos(th) + x;
yunit4 = r * sin(th) + y;
h = plot(xunit4, yunit4,'k');
hold off
end
function [AMatch,ThetaMatch,Gammas,GammaSorLoptions] = MatchingCircle(x,y,r)
% [AMatch,ThetaMatch,Gammas,GammaSorLoptions] = MatchingCircle(real(C),imag(C),R)
hold on
th = 0:pi/50:2*pi;
xunit1 = 1 * cos(th) + 0;
yunit1 = 1 * sin(th) + 0;
h = plot(xunit1, yunit1,'r');
xunit2 = 0.5 * cos(th) + 0.5;
yunit2 = 0.5 * sin(th) + 0;
h = plot(xunit2, yunit2,'r');
xunit3 = 0.5 * cos(th) - 0.5;
yunit3 = 0.5 * sin(th) + 0;
h = plot(xunit3, yunit3,'b');
xunit4 = r * cos(th) + x;
yunit4 = r * sin(th) + y;
h = plot(xunit4, yunit4,'k');
hold off

%matching network
[xout1,yout1] = circcirc(0.5,0,0.5,x,y,r);
[xout2,yout2] = circcirc(-0.5,0,0.5,x,y,r);

[A1,Theta1] = Rect2Pold(xout1(1,1),yout1(1,1));
[A2,Theta2] = Rect2Pold(xout1(1,2),yout1(1,2));
[A3,Theta3] = Rect2Pold(xout2(1,1),yout2(1,1));
[A4,Theta4] = Rect2Pold(xout2(1,2),yout2(1,2));

%VSWR
[A5,Theta5] = Rect2Pold(r*cos(0)+x,r*sin(0)+y);
[A6,Theta6] = Rect2Pold(r*cos(pi/2)+x,r*sin(pi/2)+y);
[A7,Theta7] = Rect2Pold(r*cos(pi)+x,r*sin(pi)+y);
[A8,Theta8] = Rect2Pold(r*cos(3*pi/2)+x,r*sin(3*pi/2)+y);
GammaSorL1=Pol2Rectd(A5,Theta5);
GammaSorL2=Pol2Rectd(A6,Theta6);
GammaSorL3=Pol2Rectd(A7,Theta7);
GammaSorL4=Pol2Rectd(A8,Theta8);
GammaSorLoptions=[GammaSorL1;GammaSorL2;GammaSorL3;GammaSorL4];

Gammas = [A1 Theta1;A2 Theta2;A3 Theta3;A4 Theta4;A5 Theta5;A6 Theta6;A7 Theta7;A8 Theta8];
Amin=min(Gammas(1:4,1));
if A1==Amin
    AMatch = A1;
    ThetaMatch = Theta1;
elseif A2==Amin
    AMatch = A2;
    ThetaMatch = Theta2;
elseif A3==Amin
    AMatch = A3;
    ThetaMatch = Theta3;
elseif A4==Amin
    AMatch = A4;
    ThetaMatch = Theta4;
end
    
fprintf("\nMatching Network Options:\nGamma option 1:\t%f @ %f\nGamma option 2:\t%f @ %f\nGamma option 3:\t%f @ %f\nGamma option 4:\t%f @ %f\n",A1,Theta1,A2,Theta2,A3,Theta3,A4,Theta4);
fprintf("\nGammaSorL Options for VSWR:\nGamma option 1:\t%f @ %f\nGamma option 2:\t%f @ %f\nGamma option 3:\t%f @ %f\nGamma option 4:\t%f @ %f\n",A5,Theta5,A6,Theta6,A7,Theta7,A8,Theta8);

end
function [C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut)
% S11=Pol2Rectd(0.56,-78);
% S12=Pol2Rectd(0.05,33); 
% S21=Pol2Rectd(8.64,122);
% S22=Pol2Rectd(0.66,-42);
% InOut="In"; %"In" or "Out"
% [C, R] = StabilityCircleCalc(S11,S12,S21,S22,InOut);

Del = S11*S22-S12*S21;
if InOut=="Out"
    C = conj(S22-Del*conj(S11))/(norm(S22)^2-norm(Del)^2);
    R = norm(S12*S21/(norm(S22)^2-norm(Del)^2));
elseif InOut=="In"
    C = conj(S11-Del*conj(S22))/(norm(S11)^2-norm(Del)^2);
    R = norm(S12*S21/(norm(S11)^2-norm(Del)^2));
end

[A,Theta] = Rect2Pold(real(C),imag(C));
fprintf("\nCenter = %f @ %f\tor\t%f + j%f\nRadius = %f\n",A,Theta,real(C),imag(C),R)
end
function [S11,S12,S21,S22,S]=T2S(T11,T12,T21,T22)
% T11=Pol2Rectd(0.75,-160);
% T12=Pol2Rectd(0,0);
% T21=Pol2Rectd(5,90);
% T22=Pol2Rectd(0.2,-30);
% [S11,S12,S21,S22,S]=T2S(T11,T12,T21,T22)
S11 = T21/T11;
S12 = T22 - T21*T12/T11;
S21 = 1/T11;
S22 = -T12/T11;
S = [S11 S12;
    S21 S22];
end
function [T11,T12,T21,T22,T]=S2T(S11,S12,S21,S22)
% S11=Pol2Rectd(0.75,-160);
% S12=Pol2Rectd(0,0); 
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% [T11,T12,T21,T22,T]=S2T(S11,S12,S21,S22)
T11 = 1/S21;
T12 = -S22/S21;
T21 = S11/S21;
T22 = S12 - S11*S22/S21;
T = [T11 T12;
    T21 T22];
end
function [Gammai,Gilin,gi, C, R] = GainCircleCalc(GidB,Sii)
% GidB=2; 
% Sii=Pol2Rectd(0.5,-90);
% [Gammai,Gilin,gi, C, R] = GainCircleCalc(GidB,Sii);

[A,d]=Rect2Pold(real(Sii),imag(Sii));
Gilin=10^(GidB/10);
gi=Gilin*(1-norm(Sii^2));

%| GammaS - Center | = Radius or | GammaL - Center | = Radius 
C = (gi*conj(Sii))/(1-(norm(Sii)^2)*(1-gi));
[MagC,ThetaC] = Rect2Pold(real(C),imag(C));
% ThetaCd=ThetaC*180/pi;
R = (sqrt(1-gi)*(1-norm(Sii)^2))/(1-(norm(Sii)^2)*(1-gi));

syms Gamma

f = Gamma-C==R;
Gammai(1,1)=solve(f);
[MagGammai1,ThetaGammai1] = Rect2Pold(real(Gammai(1,1)),imag(Gammai(1,1)));
% ThetaGammai1d=ThetaGammai1*180/pi;

f = Gamma-C==-R;
Gammai(1,2)=solve(f);
[MagGammai2,ThetaGammai2] = Rect2Pold(real(Gammai(1,2)),imag(Gammai(1,2)));
% ThetaGammai2d=ThetaGammai2*180/pi;

fprintf("\nCenter = %f + j%f or %fe^j%f",real(C),imag(C),MagC,ThetaC);
fprintf("\nRadius = %f\n",R);
% fprintf("\nGammai = %f + j%f or %fe^j%f",real(Gammai(1,1)),imag(Gammai(1,1)),MagGammai1,ThetaGammai1);
% fprintf("\nGammai = %f + j%f or %fe^j%f\n",real(Gammai(1,2)),imag(Gammai(1,2)),MagGammai2,ThetaGammai2);

end
function [gi] = NormGainFactCalc(Gammai,Sii)
% Sii=Pol2Rectd(0.5,-90); %Max condition |Sii| < 1 & can be either S11 or S22
% Gammai=conj(Sii); %Max condition when Gammai=conj(Sii) & can be either GammaL or GammaS
% [gi] = NormGainFactCalc(Gammai,Sii)

gi=(1-norm(Gammai)^2)*(1-norm(Sii)^2)/norm(1-Sii*Gammai)^2;
gidB = 10*log10(gi);
fprintf("\ngi = %f (linear) or %f (dB)\n",gi,gidB);

end
function [MS, ML, PAVS, PIN, PL, PAVN] = PowerCalc(GammaS,GammaL,GammaIN,GammaOUT,E1,Z1,GTlin)
% GammaS = Pol2Rectd(2.5,30);
% GammaL = Pol2Rectd(2.5,30);
% GammaIN = Pol2Rectd(2.5,30);
% GammaOUT = Pol2Rectd(2.5,30);
% E1 = 10;  %Voltage
% Z1 = 50;
% GT = 9.5;
% [MS, ML, PAVS, PIN, PL, PAVN] = PowerCalc(GammaS,GammaL,GammaIN,GammaOUT,E1,Z1,GTlin);

MS = (1-norm(GammaS)^2)*(1-norm(GammaIN)^2)/norm(1-GammaIN*GammaS)^2;
MSdB = 10*log10(MS);
ML = (1-norm(GammaL)^2)*(1-norm(GammaOUT)^2)/norm(1-GammaOUT*GammaL)^2;
MLdB = 10*log10(ML);

PAVS = norm(E1)^2/(8*real(Z1));
PIN = PAVS*MS;
PL = GTlin*PAVS;
PAVN = PL/ML;

fprintf("\nMS = %f\t(Linear)\nMS = %f\t(dB)\nML = %f\t(Linear)\nML = %f\t(dB)\n\nPAVS = %f W\nPIN = %f W\nPL = %f W\nPAVN = %f W\n",MS,MSdB,ML,MLdB,PAVS,PIN,PL,PAVN);

MagGammaa = sqrt(1-MS);
MagGammab = sqrt(1-ML);

VSWRIN = (1+MagGammaa)/(1-MagGammaa);
VSWROUT = (1+MagGammab)/(1-MagGammab);

fprintf("\n|Gamma a| = %f\n|Gamma b| = %f\nVSWRIN = %f\nVSWROUT = %f\n",MagGammaa,MagGammab,VSWRIN,VSWROUT)

end
function [GAlin, GAdB, G0, GS] = UniGACalc(S21,S11,GammaOUT,GammaS)
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% GammaOUT=Pol2Rectd(0.2,-30);
% GammaS=S11; %Max condition S22
% [GAlin, GAdB, G0, GS] = UniGACalc(S21,S11,GammaOUT,GammaS);

[GS] = UniGSCalc(GammaS,S11);
[G0] = UniG0Calc(S21);
[GL] = 1/(1-norm(GammaOUT)^2);

GSdB = 10*log10(GS);
G0dB = 10*log10(G0);
GLdB = 10*log10(GL);

GAlin = GS*G0*GL;
GAdB = 10*log10(GAlin);

fprintf("\nGS = %f (linear)\tor\t%f (dB)\nG0 = %f (linear)\tor\t%f (dB)\nGL = %f (linear)\tor\t%f (dB)\n\nGA = %f\t(Linear)\tor\t%f\t(dB)\n",GS,GSdB,G0,G0dB,GL,GLdB,GAlin,GAdB);

end
function [GPlin, GPdB, G0, GS] = UniGPCalc(S21,S22,GammaIN,GammaL)
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% GammaIN=Pol2Rectd(0.2,-30);
% GammaL=S22; %Max condition S22
% [GPlin, GPdB, G0, GS] = UniGPCalc(S21,S22,GammaIN,GammaL);

[GS] = 1/(1-norm(GammaIN)^2);
[G0] = UniG0Calc(S21);
[GL] = UniGLCalc(GammaL,S22);

GSdB = 10*log10(GS);
G0dB = 10*log10(G0);
GLdB = 10*log10(GL);

GPlin = GS*G0*GL;
GPdB = 10*log10(GPlin);

fprintf("\nGS = %f (linear)\tor\t%f (dB)\nG0 = %f (linear)\tor\t%f (dB)\nGL = %f (linear)\tor\t%f (dB)\n\nGP = %f\t(Linear)\tor\t%f\t(dB)\n",GS,GSdB,G0,G0dB,GL,GLdB,GPlin,GPdB);

end
function [GTlin, GTdB, GammaIN, GammaOUT, GL, G0, GS] = UniGTCalc(S11,S12,S21,S22,GammaS,GammaL,type)
% S11=Pol2Rectd(0.75,-160);
% S12=Pol2Rectd(0,0); %Unilateral when S12=0
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% GammaS=conj(S11); %Max condition conj(S11)
% GammaL=conj(S22); %Max condition conj(S22)
% type="UNI";  %"IN","OUT","UNI"
% [GTlin, GTdB, GammaIN, GammaOUT, GL, G0, GS] = UniGTCalc(S11,S12,S21,S22,GammaS,GammaL,type);

if norm(S11) >= 1
    fprintf("\n********************************************************************\nWARNING:\tS11 IS MAKING THE SYSTEM CONDITIONALLY STABLE\n********************************************************************\n");
end
if norm(S22) >= 1
    fprintf("\n********************************************************************\nWARNING:\tS22 IS MAKING THE SYSTEM CONDITIONALLY STABLE\n********************************************************************\n");
end

[Del K]=Del_K(S11,S12,S21,S22);

[GammaIN] = UniGammaINCalc(S11,S12,S21,S22,GammaL);
[MagGammaIN,ThetaGammaIN] = Rect2Pold(real(GammaIN),imag(GammaIN));
% ThetadGammaIN=ThetaGammaIN*180/pi;
[GammaOUT] = UniGammaOUTCalc(S11,S12,S21,S22,GammaS);
[MagGammaOUT,ThetaGammaOUT] = Rect2Pold(real(GammaOUT),imag(GammaOUT));
% ThetadGammaOUT=ThetaGammaOUT*180/pi;

if type == "IN"
    [GS] = UniGSCalc(GammaS,GammaIN);
    [G0] = UniG0Calc(S21);
    [GL] = UniGLCalc(GammaL,S22);
elseif type == "OUT"
    [GS] = UniGSCalc(GammaS,S11);
    [G0] = UniG0Calc(S21);
    [GL] = UniGLCalc(GammaL,GammaOUT);
elseif type == "UNI"
    [GS] = UniGSCalc(GammaS,S11);
    [G0] = UniG0Calc(S21);
    [GL] = UniGLCalc(GammaL,S22);
end

GSdB = 10*log10(GS);
G0dB = 10*log10(G0);
GLdB = 10*log10(GL);

GTlin = GS*G0*GL;
GTdB = 10*log10(GTlin);

fprintf("\nGS = %f (linear)\tor\t%f (dB)\nG0 = %f (linear)\tor\t%f (dB)\nGL = %f (linear)\tor\t%f (dB)\n\nGammaIN = %f + j%f = %fe^j%f\nGammaOUT = %f + j%f = %fe^j%f\n\nGT = %f + j%f\t(Linear)\tor\tGT = %f + j%f\t(dB)\n",GS,GSdB,G0,G0dB,GL,GLdB,real(GammaIN),imag(GammaIN),MagGammaIN,ThetaGammaIN,real(GammaOUT),imag(GammaOUT),MagGammaOUT,ThetaGammaOUT,real(GTlin),imag(GTlin),real(GTdB),imag(GTdB));

end
function [GammaIN] = UniGammaINCalc(S11,S12,S21,S22,GammaL)
% S11=Pol2Rectd(0.75,-160);
% S12=Pol2Rectd(0,0); %Unilateral when S12=0
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% GammaL=S22; %Max condition S22
% [GammaIN] = UniGammaINCalc(S11,S12,S21,S22,GammaL);

GammaIN=S11+(S12*S21*GammaL)/(1-S22*GammaL);

end
function [GammaOUT] = UniGammaOUTCalc(S11,S12,S21,S22,GammaS)
% S11=Pol2Rectd(0.75,-160);
% S12=Pol2Rectd(0,0); %Unilateral when S12=0
% S21=Pol2Rectd(5,90);
% S22=Pol2Rectd(0.2,-30);
% GammaS=S11; %Max condition S11
% [GammaOUT] = UniGammaOUTCalc(SS11,S12,S21,S22,GammaS);

GammaOUT=S22+(S12*S21*GammaS)/(1-S11*GammaS);

end
function [GL] = UniGLCalc(GammaL,GammaOUT)
% S22 = Pol2Rectd(0.75,-160);
% GammaL = Pol2Rectd(0.75,-160);
% [GL] = UniGLCalc(GammaL,GammaOUT);

GL=(1-norm(GammaL)^2)/(norm(1-GammaOUT*GammaL)^2);

end
function [G0] = UniG0Calc(S21)
% S21 = Pol2Rectd(0.75,-160);
% [G0] = UniG0Calc(S21)

G0=norm(S21)^2;

end
function [GS] = UniGSCalc(GammaS,GammaIN)
% S11 = Pol2Rectd(0.75,-160);
% GammaS = Pol2Rectd(0.75,-160);
% [Gs] = UniGSCalc(GammaS,GammaIN)

GS=(1-norm(GammaS)^2)/(norm(1-GammaIN*GammaS)^2);

end
function [Radius, Center, C_Angle]=Unilateral_Stability(S11_Cart, S12_Cart, S21_Cart, S22_Cart,Gs,GL)
gS = Gs*(1-abs(S11_Cart)^2)
gL = GL*(1-abs(S22_Cart)^2)
disp("L = Output S = Input")
if Circle_Wanted == "Input"
    disp("Input Circle")
    B_S = 1-(abs(S11_Cart(1))^2*(1-gS));
    Radis = (sqrt(1-gS)*(1-abs(S11_Cart)^2))/(B_S)
    Center = (gS*(S11_Cart(1)))/(B_S)
    C_Angle = -S11_Cart(2)
elseif Circle_Wanted == "Output"
    disp("Output Circle")
    B_L = 1-(abs(S22_Cart(1))^2*(1-gL));
    Radius = (sqrt(1-gL)*(1-abs(S22_Cart)^2))/(B_L)
    Center = (gL*(S22_Cart(1)))/(B_L)
    C_Angle = -S22_Cart(2)
end
end
function [G_sMax,G_LMax, G_0, GTU_MaxdB] = TransducerGain(S11_Cart,S12_Cart,S21_Cart,S22_Cart)
G_sMax = 1/(1-norm(S11_Cart)^2);
G_LMax = 1/(1-norm(S22_Cart)^2);
G_0 = norm(S21_Cart)^2;
GTU_Max = G_sMax*G_0*G_LMax;
GTU_MaxdB = 10*log10(GTU_Max);
end
function StabilityCheck(S11,S12,S21,S22)

%inputs
% S11 = Pol2Rectd(0.75,-160);
% S12 = Pol2Rectd(0,0);
% S21 = Pol2Rectd(5,90);
% S22 = Pol2Rectd(0.2,-30);
% fc=3e9;
% StabilityCheck(S11,S12,S21,S22)

%functions
[Del, K]=Del_K(S11,S12,S21,S22);
[RadiusL, CenterL, RadiusS, CenterS]=LSRadCent(S11,S12,S21,S22,Del);
[GTU_Max, GTU_MaxdB]=MaxGain(S11,S21,S22);

end
function [Del, K]=Del_K(S11,S12,S21,S22)
%Inputs
% S11 = Pol2Rectd(0.75,-160);
% S12 = Pol2Rectd(0,0);
% S21 = Pol2Rectd(5,90);
% S22 = Pol2Rectd(0.2,-30);
% [Del K]=Del_K(S11,S12,S21,S22)

%Functions
Del = (S11)*(S22) - (S12)*(S21);
[A,Theta] = Rect2Pold(real(Del),imag(Del));
K = (1-(norm(S11)^2)-(norm(S22)^2)+norm(Del)^2)/(2*norm((S12)*(S21)));
fprintf("\nDel = %f @ %f\tor\t%f + j%f\n|Del| = %f\nK = %f\n",A,Theta,real(Del),imag(Del),norm(Del),K);
if K > 1 && norm(Del) < 1
    fprintf("\nThis is unconditionally stable\n")
else
    fprintf("\nThis is potentially unstable\nThe origin on the GammaL plane (Input) is stable if S11 = %f < 1\nThe origin on the GammaL plane (Output) is stable if S22 = %f < 1\n",norm(S11),norm(S22));
end

end
function [RadiusL, CenterL, RadiusS, CenterS]=LSRadCent(S11,S12,S21,S22,Del)
%Inputs
% S11 = Pol2Rectd(0.75,-160);
% S12 = Pol2Rectd(0,0);
% S21 = Pol2Rectd(5,90);
% S22 = Pol2Rectd(0.2,-30);
% Del = 2;
% [RadiusL CenterL RadiusS CenterS]=LSRadCent(S11,S12,S21,S22,Del);

%Functions
RadiusL = norm((S12)*(S21)/(norm(S22)^2-norm(Del)^2));
CenterL = conj((S22)-(Del*conj(S11)))/(norm(S22)^2-norm(Del)^2);
RadiusS = norm((S12)*(S21)/(norm(S11)^2-norm(Del)^2));
CenterS = conj((S11)-(Del*conj(S22)))/(norm(S11)^2-norm(Del)^2);
fprintf("\nL = Output S = Input\nRadiusL = %f\t&\tCenterL = %f + j%f\nRadiusS = %f\t&\tCenterS = %f + j%f\n",RadiusL,real(CenterL),imag(CenterL),RadiusS,real(CenterS),imag(CenterS));

end
function [GTU_Max, GTU_MaxdB]=MaxGain(S11,S21,S22)
%Inputs
% S11 = Pol2Rectd(0.75,-160);
% S21 = Pol2Rectd(5,90);
% S22 = Pol2Rectd(0.2,-30);
% [GTU_Max GTU_MaxdB]=MaxGain(S11,S21,S22)

%Functions
G_sMax = 1/(1-norm(S11)^2);
G_LMax = 1/(1-norm(S22)^2);
G_0 = norm(S21)^2;
GTU_Max = G_sMax*G_0*G_LMax;
GTU_MaxdB = 10*log10(GTU_Max);
fprintf("\nGTU_Max = %f\t&\tGTU_MaxdB = %f\n",GTU_Max,GTU_MaxdB);

end
function EQNFinder(Detail)
%This function looks up eqations based on the detail string input
%Detail = "GTU_MAX"; 
%EQNFinder(Detail);
syms ML MS GP PIN PAVN PL PAVS GA GT RadiusS CenterS CenterL RadiusL GTU_Max GammaL GammaIN GammaOUT GammaS GammaIN S11 S12 S21 S22 Delta K G_sMax G_LMax G_0

Detail=upper(Detail);

if Detail == "CENTERS" || Detail == "ALL"
    pretty(CenterS == conj((S11)-(Delta*(S22)))/(norm(S11)^2-norm(Delta)^2));
end
if Detail == "RADIUSS" || Detail == "ALL"
    pretty(RadiusS == abs((S12)*(S21)/(norm(S11)^2-norm(Delta)^2)));
end
if Detail == "CENTERL" || Detail == "ALL"
    pretty(CenterL == conj((S22)-(Delta*(S11)))/(norm(S22)^2-norm(Delta)^2));
end
if Detail == "RADIUSL" || Detail == "ALL"
    pretty(RadiusL == abs((S12)*(S21)/(abs(S22)^2-abs(Delta)^2)));
end
if Detail == "GTU_MAX" || Detail == "ALL"
    pretty(G_sMax == 1/(1-norm(S11)^2));
    pretty(G_LMax == 1/(1-norm(S22)^2));
    pretty(G_0 == norm(S21)^2);
    pretty(GTU_Max == 1/(1-norm(S11)^2)*norm(S21)^2*1/(1-norm(S22)^2));
end
if Detail == "K" || Detail == "ALL"
    pretty(K==(1-norm(S11)^2-norm(S22)^2+norm(Delta)^2)/(2*abs(S12*S21)));
    fprintf("If K > 1 then the amp is unconditionally stable\n\n");
end
if Detail == "DELTA" || Detail == "ALL"
    pretty(Delta==S11*S22-S12*S21);
    fprintf("If |Delta| < 1 then the amp is unconditionally stable\n\n");
end
if Detail == "GAMMAIN" || Detail == "ALL"
    pretty(GammaIN == S11+(S12*S21*GammaL)/(1-S22*GammaL));
end
if Detail == "GAMMAOUT" || Detail == "ALL"
    pretty(GammaOUT == S22+(S12*S21*GammaS)/(1-S11*GammaS));
end
if Detail == "GT" || Detail == "GTU" || Detail == "ALL"
    pretty(GT == (1-norm(GammaS)^2)/(norm(1-S11*GammaS)^2)*norm(S21)^2*(1-norm(GammaL)^2)/(norm(1-S22*GammaL)^2));
    pretty(GT == PL/PAVS);
end
if Detail == "GA" || Detail == "ALL"
    pretty(GA == (1-norm(GammaS)^2)/(norm(1-S11*GammaS)^2)*norm(S21)^2*1/(1-norm(GammaOUT)^2));
    pretty(GA == PAVN/PAVS);
end
if Detail == "GP" || Detail == "ALL"
    pretty(GP == (1/(1-norm(GammaIN)^2))*norm(S21)^2*(1-norm(GammaL)^2)/(norm(1-S22*GammaL)^2));
    pretty(GP == PL/PIN);
end
if Detail == "MS" || Detail == "ALL"
    pretty(MS == (1-norm(GammaS)^2)*(1-norm(GammaIN)^2)/norm(1-GammaIN*GammaS)^2);
end
if Detail == "ML" || Detail == "ALL"
    pretty(ML == (1-norm(GammaL)^2)*(1-norm(GammaOUT)^2)/norm(1-GammaOUT*GammaL)^2);
end


end
function RF2_HW04_Prob4()
%Inputs
S11 = Pol2Rectd(0.3,-70);
S21 = Pol2Rectd(3.5,85);
S12 = Pol2Rectd(0.2,-10);
S22 = Pol2Rectd(0.4,-45);
Vs=5;
Zs=40;
Zl=73;
Z0=50;

%Formula
%a)
GammaS=(Zs-Z0)/(Zs+Z0);
GammaL=(Zl-Z0)/(Zl+Z0);
GammaIN=S11+(S12*S21*GammaL)/(1-S22*GammaL);
GT=((1-abs(GammaS)^2)/abs(1-GammaIN*GammaS)^2)*(abs(S21)^2)*((1-abs(GammaL)^2)/abs(1-S22*GammaL)^2);
GTdB=10*log10(GT);
fprintf("\na)\nGT = %f (or %f dB)\n",GT,GTdB);
%b)
GP=((1)/(1-abs(GammaIN)^2))*(abs(S21)^2)*((1-abs(GammaL)^2)/abs(1-S22*GammaL)^2);
GPdB=10*log10(GP);
fprintf("\nb)\nGP = %f (or %f dB)\n",GP,GPdB);
%c)
GammaOUT=S22+(S12*S21*GammaS)/(1-S22*GammaS);
GA=((1-abs(GammaS)^2)/abs(1-GammaIN*GammaS)^2)*(abs(S21)^2)*((1)/(1-abs(GammaOUT)^2)); 
GAdB=10*log10(GA);
fprintf("\nc)\nGA = %f (or %f dB)\n",GA,GAdB);
%d)
PAVS=(Vs^2)/(8*real(Zs));
PAVSdB=10*log10(PAVS);
fprintf("\nd)\nPAVS = %f W (or %f dBW)\n",PAVS,PAVSdB);
%e)
PINC=PAVS*(1-abs(S11)^2);
PINCdB=10*log10(PINC);
fprintf("\ne)\nPINC = %f W (or %f dBW)\n",PINC,PINCdB);
%f)
PREF=PAVS-PINC;
PREFdB=10*log10(PREF);
fprintf("\nf)\nPREF = %f W (or %f dBW)\n",PREF,PREFdB);
%g)
PIN=PAVS*(1-abs(GammaS)^2)*(1-abs(GammaIN)^2)/abs(1-GammaS*GammaIN)^2;
PINdB=10*log10(PIN);
fprintf("\ng)\nPIN = %f W (or %f dBW)\n",PIN,PINdB);
%h)
PAVN=GA*PAVS;
PAVNdB=10*log10(PAVN);
fprintf("\nh)\nPAVN = %f W (or %f dBW)\n",PAVN,PAVNdB);
%i)
PL=PAVN*(1-abs(GammaL)^2)*(1-abs(GammaOUT)^2)/abs(1-GammaOUT*GammaL)^2;
PLdB=10*log10(PL);
fprintf("\ni)\nPL = %f W (or %f dBW)\n",PL,PLdB);

fprintf("\n\t***************************\n");
fprintf("\n\t*****  Formulas used  *****\n");
fprintf("\n\t***************************\n");

syms GammaS GammaIN GammaS S11 S12 S21 GammaOUT S22 GammaL GT GP GA PAVS PINC PREF PIN PAVN PL Zs Vs
pretty(GT==((1-abs(GammaS)^2)/abs(1-GammaIN*GammaS)^2)*(abs(S21)^2)*((1-abs(GammaL)^2)/abs(1-S22*GammaL)^2));
pretty(GP==((1)/(1-abs(GammaIN)^2))*(abs(S21)^2)*((1-abs(GammaL)^2)/abs(1-S22*GammaL)^2));
pretty(GA==((1-abs(GammaS)^2)/abs(1-GammaIN*GammaS)^2)*(abs(S21)^2)*((1)/(1-abs(GammaOUT)^2)));
pretty(PAVS==(Vs^2)/(8*real(Zs)));
pretty(PINC==PAVS*(1-abs(S11)^2));
pretty(PREF==PAVS-PINC);
pretty(PIN==PAVS*(1-abs(GammaS)^2)*(1-abs(GammaIN)^2)/abs(1-GammaS*GammaIN)^2);
pretty(PAVN==GA*PAVS);
pretty(PL==PAVN*(1-abs(GammaL)^2)*(1-abs(GammaOUT)^2)/abs(1-GammaOUT*GammaL)^2);
end
function RF2_HW02_Prob3()
% RF2_HW02_Prob3()
syms L C a3 b3 a4 b4

% disp("Problem 1)");
% disp(" ");
% 
% f1=0.97+1i*0.21==a3+b3;
% 
% f2=0.89-1i*0.41==a3-b3;
% 
% [a3 b3]=solve(f1,f2)

disp("Problem 3)");
disp(" ");
f=500e6;
f1=1i*50==(1i*2*pi*f*L * 1/(1i*2*pi*f*C))/(1i*2*pi*f*L + 1/(1i*2*pi*f*C));

f=2e9;
f2=-1i*50==(1i*2*pi*f*L * 1/(1i*2*pi*f*C))/(1i*2*pi*f*L + 1/(1i*2*pi*f*C));

[C, L]=solve(f1,f2);

disp("C is ");
disp(C);
disp("L is ");
disp(L);

% disp("Problem 4)");
% disp(" ");
% 
% Z0=50;
% f=2e9;
% w=2*pi*f;
% 
% c2=1.59e-12;
% z2=1/(1i*w*c2*Z0)
% 
% l3=3.98e-9;
% Z3=1i*w*l3/Z0
% 
% l4=3.98e-9;
% Z4=1i*w*l4/Z0
% 
% c5=1.68e-12;
% z5=1/(1i*w*c5*Z0)
% 
% l6=1.99e-9;
% Z6=1i*w*l6/Z0

end
function Extra()

Z0=50;
R1=5;
R2=45;
L=1e-9;
C=1e-12;

fr=1/(2*pi*sqrt(L*C));

% f=5.0329e9;
f=[0:100e6:100e9];
w=2*pi.*f;

ZL=1j.*w.*L;
ZC=1./(1j.*w.*C);

Z1=R2.*ZC./(R2+ZC);
Z2=Z1.*ZL./(Z1+ZL);
Zin=Z2+R1;
Zin=Zin/Z0;
disp(Zin);
Zin_real=real(Zin);
Zin_imag=imag(Zin);
plot(f,Zin_real,f,Zin_imag)

end
function RF2_HW01_Prob1()
syms lambda
n=[1:0.05:32];
z=lambda./n;

Z0=50;
a=0.3;
b=0.4;
beta=2*pi/lambda;
gammaL=a+b*1j;
rho=abs(gammaL);
theta=angle(gammaL);
gammaL=rho*exp(1j*theta);
gammaL_betapos=gammaL
Zin_betapos=Z0*(1+gammaL)/(1-gammaL)
V=gammaL*exp(1j*(-1)*beta*z)+exp(1j*beta*z);
Vmag1=round(abs(V),4);
I=(gammaL*exp(1j*(-1)*beta*z)-exp(1j*beta*z));
Imag1=round(abs(I),4);
VImag1=Vmag1./Imag1;

b=-0.4;
gammaL=a+b*1j;
rho=abs(gammaL);
theta=angle(gammaL);
gammaL=rho*exp(1j*theta);
gammaL_betapneg=gammaL
Zin_betaneg=Z0*(1+gammaL)/(1-gammaL)
V=gammaL*exp(1j*(-1)*beta*z)+exp(1j*beta*z);
Vmag2=round(abs(V),4);
I=(gammaL*exp(1j*(-1)*beta*z)-exp(1j*beta*z));
Imag2=round(abs(I),4);
VImag2=Vmag2./Imag2;

figure(1)
subplot(3,2,1);
plot(beta*z,Vmag1);subtitle('|V| with +b');xlabel('beta*z');
subplot(3,2,3);
plot(beta*z,Imag1);subtitle('|I| with +b');xlabel('beta*z');
subplot(3,2,5);
plot(beta*z,VImag1);subtitle('|V|/|I| with +b');xlabel('beta*z');
subplot(3,2,2);
plot(beta*z,Vmag2);subtitle('|V| with -b');xlabel('beta*z');
subplot(3,2,4);
plot(beta*z,Imag2);subtitle('|I| with -b');xlabel('beta*z');
subplot(3,2,6);
plot(beta*z,VImag2);subtitle('|V|/|I| with -b');xlabel('beta*z');
end
function RF2_HW01_Prob2()
%Not working
syms lambda s11 s12 s21 s22 a1 a2 a12 b1 b2 b12
beta=2*pi/lambda;

Z0=50;
Z0_left=Z0;
Z0_right=Z0;
Z01=50;
Z02=100;
ZL=200;
Length01=lambda/8;
Length02=lambda/4;

disp("*** Find S11 ***");
disp(" ");
Z1=(Z02^2)/Z0_right;
Z2=ZL+Z1;
Zin=Z0_left*(Z2+1j*Z0_left*tan(beta*Length01))/(Z0_left+1j*Z2*tan(beta*Length01));
disp("Zin=");
disp(round(Zin,4));
disp("S11=|Gamma|=|(Zin-Z0)/(Zin+Z0)|=");
Gamma=(Zin-Z0_left)/(Zin+Z0_left);
disp(round(abs(Gamma),4));

disp("*** Find S21 ***");
% disp(" ");

% a1=1/sqrt(Z0)
% 
% mag2=sqrt(Z01);
% a2=a1*exp(1j*(-1)*beta*Length01);
% b2=b1*exp(1j*beta*Length01);
% v1=mag2*(a2 + b2);
% magv2=mag2*(Z1/(ZL+Z1));
% v2=magv2*v1;
% mag3=magv2;
% a3=a2*exp(1j*(-1)*beta*Length02);
% b3=b2*exp(1j*beta*Length02);
% mag4=mag3*sqrt(Z02);
% a4=a3;
% b4=b3;
% round(mag4*(a4+b4),4);
% [a2 b2] = solve(mag4*(a4+b4)==-sqrt(Z0)*(a2+b2),(1/sqrt(Z02))*mag4*(a4+b4)==-(1/sqrt(Z0))*(a2+b2))

a1=1/sqrt(50);
b1=0.7778/sqrt(50);
% [a2 b2] = solve(5*a1*exp(1j*(-1)*3*pi/4)+5*b1*exp(1j*3*pi/4)==(a2+b2),2.5*a1*exp(1j*(-1)*3*pi/4)-2.5*b1*exp(1j*3*pi/4)==(b2-a2));
% mag_a2=round(a2,4);
% mag_b2=round(b2,4);
% mag_a2=round(abs(a2),4);
% mag_b2=round(abs(b2),4);

% a1=1;
% b1=1;
f1=(sqrt(50)/2)*(a1*exp(1j*(-1)*pi/4)+b1*exp(1j*pi/4))==sqrt(100)*(a12+b12);
f2=(1/sqrt(50))*(a1*exp(1j*(-1)*pi/4)-b1*exp(1j*pi/4))==(1/sqrt(100))*(a12-b12);
[a12, b12]=solve(f1,f2);
% a12round=round(a12,4)
% b12round=round(b12,4)

a13=a12*exp(1j*(-1)*pi/2);
b13=b12*exp(1j*pi/2);
% a13round=round(a13,4)
% b13round=round(b13,4)


f1=sqrt(100)*(  (5/8)*sqrt(2)*a1*exp(1j*(-1)*3*pi/4)    -   (3/8)*sqrt(2)*b1*exp(1j*(-1)*pi/4)  -   (3/8)*sqrt(2)*a1*exp(1j*pi/4)   +   (5/8)*sqrt(2)*b1*exp(1j*3*pi/4) ==  sqrt(50)*(0+b2));
f2=(1/sqrt(100))*(  (5/8)*sqrt(2)*a1*exp(1j*(-1)*3*pi/4)  -   (3/8)*sqrt(2)*b1*exp(1j*(-1)*pi/4)    +   (3/8)*sqrt(2)*a1*exp(1j*pi/4)   -   (5/8)*sqrt(2)*b1*exp(1j*3*pi/4) ==  -(1/sqrt(50))*(0-b2));
% [a2 b2]=solve(f1,f2);
[b2]=solve(f2);
% a2round=round(a2,4)
b2round=round(b2,4);
% a2round=round(abs(a2),4)
b2round=round(abs(b2),4);

% f1=-0.089-1j*0.011;
% f2=;
% [a12 b12]=solve(f1,f2);

% a1=1/sqrt(50);
% b1=0.7778/sqrt(50);
% a2=abs(-0.42+1i*0.17);
% b2=abs(-0.47-1i*0.28);
% a1=1/sqrt(50);
% b1=0.7778/sqrt(50);
a2=abs(a2);
b2=abs(b2);
[s11, s12 s21 s22]=solve(b1==s11*a1+s12*a2,b2==s21*a1+s22*a2);
s11round = round(s11,4)
s12round = round(s12,4);
s21round = round(s21,4)
s22round = round(s22,4);

disp("I know there is some issue with my s21 value compared to the value"); 
disp("calculated using the ABCD parameters. I'm not certain as to what the");
disp("issue is unfortunately.");
disp(" ");

%************************************************************
%Verification steps
%************************************************************

disp("*** Verification ***");
%TL01
TL01_A=cos(beta*Length01);
TL01_B=1i*Z01*sin(beta*Length01);
TL01_C=1i*(1/Z01)*sin(beta*Length01);
TL01_D=cos(beta*Length01);
TL01_ABCD=[TL01_A TL01_B;TL01_C TL01_D];
% disp("TL01_ABCD=");
% disp(round(TL01_ABCD,4));
%TLZL
TLZL_A=1;
TLZL_B=ZL;
TLZL_C=0;
TLZL_D=1;
TLZL_ABCD=[TLZL_A TLZL_B;TLZL_C TLZL_D];
% disp("TLZL_ABCD=");
% disp(round(TLZL_ABCD,4));
%TL02
TL02_A=cos(beta*Length02);
TL02_B=1i*Z02*sin(beta*Length02);
TL02_C=1i*(1/Z02)*sin(beta*Length02);
TL02_D=cos(beta*Length02);
TL02_ABCD=[TL02_A TL02_B;TL02_C TL02_D];
% disp("TL02_ABCD=");
% disp(round(TL02_ABCD,4));
%Total
Total_ABCD=TL01_ABCD*TLZL_ABCD*TL02_ABCD;
% disp("Total_ABCD=");
% disp(round(Total_ABCD,4));
%ABCD >> S
A=Total_ABCD(1,1);
B=Total_ABCD(1,2);
C=Total_ABCD(2,1);
D=Total_ABCD(2,2);
[S11,S12,S21,S22] = ABCD2S(A,B,C,D,Z0);
S=[S11,S12;S21,S22];
round(S,4);
round(abs(S),4);
disp("S11=");
disp(round(abs(S(1,1)),4));  %Mag only
disp("S21=");
disp(round(abs(S(2,1)),4));  %Mag only
end