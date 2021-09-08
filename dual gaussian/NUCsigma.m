clc;
clear;
close all;
%% calculate range in water
E=200;%MeV
p=1.735;
alpha=2.623*10^-3;
R=alpha*E^p;%cm
R=R*10;%mm
R=258;
x=1:floor(R);%mm
n=length(x);
%% calculate Wnuc
Wnuc=0.052*log(1.13+x./(11.2-0.023*R))...
    +0.35*(0.0017*R^2-R)./((R+3)^2-x.^2)...
    -1.61*10^-9.*x*(R+3)^2;

extendw=zeros(1,300-n);
Wnuc=max(0,Wnuc);
% Wnuc=[Wnuc extendw];
figure;
plot(x,Wnuc,'-')
xlabel('depth[mm]')
ylabel('Wnuc')
%% calculate sigmaNUC in[mm]
sigmaNUC=2.85+0.0014*R.*log(x+3)+0.06.*x...
    -7.4^-5.*x.^2-0.22*R./(x-R-5).^2;
figure;
plot(x,sigmaNUC,'-')
xlabel('depth[mm]')
ylabel('sigmaNUC[mm]')
