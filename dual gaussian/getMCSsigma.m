function sigmaMCS = getMCSsigma(initenergy,range)
%% calculate the Xs value
%for water give the vaule of Hydrogen and Oxygen 
% Alpha=137.3060^-1;%fine structure number
% N=6.02*10^23;%Avagadro number
% re=2.82*10^-13;%calssical electron radius[cm] 
% ZH=1;
% AH=1;
% XsH=(Alpha*N*re^2*ZH^2/AH*(2*log(33219*(AH*ZH)^(-1/3))-1))^-1;
% 
% ZO=16;
% AO=8;
% XsO=(Alpha*N*re^2*ZO^2/AO*(2*log(33219*(AO*ZO)^(-1/3))-1))^-1;
% 
% XsW=46.88*(1/9*XsH^-1+8/9*XsO^-1)^-1;
XsW=46.88;

%% calculate the fdM value
E=initenergy;%MeV
R=range;
Ep=938.272;%MeV
c=3*10^8;
tau=E/Ep;
pv1=(tau+2)/(tau+1)*E;
x=0.1:0.1:floor(R*10)/10;%cm
n=length(x);


%% calculate sigmax and sigmay
Es=15;%MeV
% TdM=fdM.*(Es./pv).^2/XsW;
pv=pv1.*sqrt(1-x./R);
% fdM=0.5244+0.1975*log10(1-(pv./pv1).^2)+0.2320.*log10(pv)-0.0098.*log10(pv).*log10(1-(pv./pv1).^2);

A=2.*(pv).^2./(R-x);
sigmax_sq=zeros(n,1);
for i=1:n
    pvhalf=pv1*sqrt(1-(x(i)/2)/R);
%     A=2*pv(i)^2/(R-x(i))
    fdM=0.5244+0.1975*log10(1-(pvhalf/pv1)^2)+0.2320*log10(pvhalf)-0.0098*log10(pvhalf)*log10(1-(pvhalf/pv1)^2);
    t=x(i)/R;
    sigmax_sq(i)=Es^2/XsW*fdM*R^2/A(i)*((1-t)^2*log((1-t)^-1)+3/2*t^2-t);
end
% t=x./R;
% sigma=(1-t).^2.*log(1./(1-t))+3/2.*t.^3-t;
% sigma=sigma';

sigmaMCS=sqrt(sigmax_sq)*10;
extend=ones(300-n,1)*sigmaMCS(n);
sigmaMCS=[sigmaMCS; extend];

end

