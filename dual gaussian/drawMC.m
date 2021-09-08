% clc;
% clear;
% close all;
m=csvread('200MeV0%2D_1.csv',8,1);
y=m(:,1);
z=m(:,2);
Dose=m(:,3);
Dose=flipud(Dose)./max(Dose);

n=length(x);
d=zeros(n,n);
for i=1:length(Dose)
    iy=y(i)+1;
    iz=z(i)+1;
    d(iz,iy)=Dose(i);
end
figure;
% d=d./max(d);
h=pcolor(d);
% axis([0 291 0 291])
dosezx=reshape(dosezx,n^2,1);
dosezx=dosezx./max(dosezx);
dosezx=reshape(dosezx,n,n);
% MCres-analyticalres
res=dosezx-d;
% max(res)
figure;
pcolor(res);
colormap jet
colorbar('location','EastOutside')
xlabel('z-axis[mm]');
ylabel('y-axis[mm]');