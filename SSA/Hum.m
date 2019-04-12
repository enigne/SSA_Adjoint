function [u,H,beta,x,fi,psi,wght,bwght,dudC]=Hum(N,xL,H0,HGL,iGL,ist,C,rhoig,a,Hcfact, m)
% analytical solution to simplified problem
% iGL is GL index, ist is index for x_*
% Hcfact is factor for the calibration at x_{GL}
%   m=1;
u=zeros(N,1);
fi=u;
psi=u;
wght=u;
bwght=u;
H=u;
x=u;
beta=u;
Np=N+1;
h=xL/N;
for i=1:N
    x(i)=i*h;
end
mm=m-1;
m1=m+1;
m2=m+2;
m3=m+3;
H0m=H0^m2;
p=1/m2;
xGL=x(iGL);
%           Hc=(m+2)*C*a^m/((m+1)*rhoig);
Hc=(H0m-HGL^m2)/xGL^m1*Hcfact;
istop=iGL;
%   istop=iGL-350;
for i=1:istop
    Hi=(H0m-Hc*x(i)^m1)^p;
    H(i)=Hi;
    u(i)=a*x(i)/Hi;
end

for i=istop+1:N
    H(i)=H(istop);
    u(i)=a*x(i)/H(i);
end
for i=1:istop
    beta(i)=C(i)*u(i)^mm;
end
xst=x(ist);
Hst=H(ist);
fic=xst/(rhoig*Hst^m3);
axstm=(a*xst)^m;
for i=1:ist
    fi(i)=-a*fic*H(i)^m;
    psi(i)=-1/Hst+C(i)*fic*((a*x(i))^m-axstm);
    wght(i)=-fi(i)*u(i)^m;
    bwght(i)=m1*C(i)*a*xst*(a*x(i))^m/(rhoig*Hst^m3*H(i));
end
dudC=a^m1*xst^m2/(m1*rhoig*Hst^m3);

