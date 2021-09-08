%% initialize
% clc;
% clear;
% close all;
E=200;%MeV
spread=.0;%percent
% alpha=2.2*10^-3;
% p=1.77;
%modifed alpha and p value add const c
alpha=2.623*10^-3;
% alpha=2.2*10^-3;
p=1.735;
% p=1.77
rho=1;
beta=0.012;
gamma=0.6;
epi=0;
%% calculate depth and sigma
R=alpha*E^p+0.25;

sigma1=0.014*R^0.935;
energyspread=E*spread;
sigma=sqrt(sigma1^2+(energyspread*alpha*p)^2*E^(2*p-2));
% sigma=0.3;
% sigma=0.02275*R+0.12085*10^-4*R^2;
% sigma=0.012*R^0.935;
%% calculate 1D dose depth curve
% x=linspace(0,10,100);
% D=zeros(1,length(x));
%proximal part
depth1=round(R-10*sigma,1);
x1=0:.1:depth1;
% y1=((R-x1).^(1/p-1)+(beta+gamma*beta*p)*(R-x1).^(1/p))./(rho*p*alpha^(1/p)*(1+beta*R));
y1=(17.93*(R-x1).^-0.435+(0.444+31.7*epi/R).*(R-x1).^0.565)/(1+0.012*R);
%distal part
depth2=round(R+5*sigma,1);
x2=depth1+0.1:.1:depth2;
y2=zeros(1,length(x2));
for i=1:length(x2)

    z=-(R-x2(i))./sigma;
    %d1=U(0.065)
    a=-1/2+0.565;
    fun1=@(t)t.^(-a-1/2).*exp(-1/2.*t.^2).*cos(z.*t+(1/2*a+1/4)*pi);
    d1=sqrt(2/pi)*exp(1/4*z^2)*integral(fun1,0,Inf);
    %d2=U(1.065)
    a=-1/2+1.565;
    fun2=@(t)t.^(a-1/2).*exp(-1/2*t.^2-z.*t);
    %gamma
    fun=@(t)t.^(1/2+a-1).*exp(-t);
    g=integral(fun,0,Inf);
    d2=exp(-1/4*z^2)/g*integral(fun2,0,Inf);
    y2(i)=(exp(-(R-x2(i))^2/(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*d1+(0.157+11/26*epi/R)*d2))/(1+0.012*R);
% y2=(exp(-(r0-z2).^2./(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*
end
%% plot 1D depth dose curve
x=zeros(1,300);
y=zeros(1,300);
len=length(x1)+length(x2);
x(1:len)=[x1 x2];
y(1:len)=[y1 y2];
figure;
y=y./max(y);
plot(x,y,'-');
xlabel('depth(cm)');
ylabel('D/Dmax');
hold on;

[maxdose,maxIdx]=max(y);
depthBP=maxIdx/10;
%% plot MC sim depth dose curve
m=csvread('200MeV0%1D_1.csv',8,2);
x1=m(:,1)./10;
y1=m(:,2);
y1=y1./max(y1);
y1=flipud(y1);
plot(x1,y1,'--','Color',[0 1 0])
hold on;
legend('analytical fit','MC sim')
% m=csvread('200MeV0%1D.csv',8,2);
% x1=m(:,1);
% y1=m(:,2);
% y1=y1./max(y1);
% x1=flipud(x1)./100;
% plot(x1,y1,'--','Color',[1 0 0])
% hold on;
%% calculate MCSsigma and NUCsigma
sigmaMCS=getMCSsigma(E,R);
[Wnuc,sigmaNUC]=getNUCsigma(E,depthBP);%depth of Bragg peak
inisigma=10;
sigmaMCS=sqrt(sigmaMCS.^2+inisigma^2);
% sigma=sqrt(sigma.^2+inisigma^2);


n=length(y);
dose=zeros(n,n);
sigma=zeros(n,1);
for i=1:n
    sigma(i)=0.002275*i+0.12085*10^-6*i^2;
end
sigma=sigma*10;
inisigma=10;
sigma=sqrt(sigma.^2+inisigma^2);
%% calculate central-axis term
ssd=2000;%mm
n=length(y);
c=zeros(n,1);
for i=1:n
    c(i)=y(i)*(ssd/(ssd+i))^2;
end
%% calculate off-axis term
centrax=n/2+.5;centray=n/2+.5;
for i=1:n
    for j=1:n
        for k=1:n
            dose(i,j,k)=c(i)*(1-Wnuc(i))*exp(-((j-centrax)^2+(k-centray)^2)/(2*sigmaMCS(i)^2))/(2*pi*sigmaMCS(i)^2)...
                +c(i)*Wnuc(i)*exp(-((j-centrax)^2+(k-centray)^2)/(2*sigmaNUC(i)^2))/(2*pi*sigmaNUC(i)^2);
            
        end
    end
end
%% 3D to 2D dose distribution
dosezx=zeros(n,n);
dosezx=sum(dose,3);
dosezx=dosezx./max(max(dosezx));
% dosezx(isnan(dosezx))=0;
figure;
h=pcolor(dosezx);











