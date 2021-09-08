function [weightNUC,sigmaNUC] = getNUCsigma(initenergy,range)
%% calculate range in water
E=initenergy;%MeV

R=range*10;%mm
x=1:floor(R);%mm
n=length(x);

%% calculate Wnuc
weightNUC=0.052*log(1.13+x./(11.2-0.023*R))...
    +0.35*(0.0017*R^2-R)./((R+3)^2-x.^2)...
    -1.61*10^-9.*x*(R+3)^2;
weightNUC=max(0,weightNUC);
extendw=zeros(1,300-n);
weightNUC=[weightNUC extendw];
%% calculate sigmaNUC in[mm]
sigmaNUC=2.85+0.0014*R.*log(x+3)+0.06.*x...
    -7.4^-5.*x.^2-0.22*R./(x-R-5).^2;
extends=ones(1,300-n)*sigmaNUC(n);
sigmaNUC=[sigmaNUC extends];
end

