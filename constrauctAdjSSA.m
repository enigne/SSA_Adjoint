function [A11, A12, A21, A22, F1, F2]=constrauctAdjSSA(N, n, ist, sigma, u, H, bxc, Afact, rhoig, dx)
% solve for the adjoint of the SSA equation
% Input
% N=number of inner points, boundaries at 0 and N+1
% n in Glen's flow law
% m in Wertman's friction law
% u velocity from forward eq
% H thickness from forward eq
% betain coeff in front of u in friction law (from forward solve)
% bxc constant inclination of bottom surface
% AA constant in viscosty (from forward solve)
% rhoig rho for ice times g
% h step length
% x0, xN first and last point of interval (not used)
% Output
% psi and fi in adjoint eq
% wght weight in perturbation integral
% A system matrix
% fic coeff multiplying fi_x in 1st eq
Cbeta=2.4126e4;
Afact=2/Afact^(1/n);
m=1/3;
% m=1;
mm1=m-1;
gamc=0;
gamc1=0;
gamc2=0;
ux=zeros(N,1);
eta=ux;
beta=ux;
fic=ux;
bx=ux;
gam=ux;
gam1=gam;
gam2=gam;
Fh=ux;
Fu=ux;
A11 = zeros(N, N);
A12 = A11;
A21 = A11;
A22 = A11;
F1 = ux;
F2 = ux;
%
p=(1-n)/n;
N1=N-1;
hh=2*dx;
h2=dx^2;
% u derivative
for i=2:N1
    ux(i)=(u(i+1)-u(i-1))/hh;
end
ux(1)=ux(2);
ux(N)=ux(N1);
% artificial damping
Ngam=ceil(0.99*N);
for i=Ngam:N
    gam(i)=gamc*(i-Ngam)/(N-Ngam);
    gam1(i)=gamc1*(i-Ngam)/(N-Ngam);
    gam2(i)=gamc2*(i-Ngam)/(N-Ngam);
end

gaussfact=1/sqrt(2*pi*sigma^2);
xi = zeros(N, 1);
for i=1:N
    xi(i)=(i-ist)*dx;
    Fu(i)=gaussfact*exp(-xi(i)^2/(2*sigma^2));
end
Fint = dx*sum(Fu);
Fu=Fu./Fint;
% viscosity, fi coeff in first eq, bottom derivative
for i=1:N
    eta(i)=Afact*abs(ux(i))^p;
    fic(i)=eta(i)*ux(i)-rhoig*H(i);
end
for i=1:ist
    beta(i)=Cbeta*abs(u(i))^mm1;
end
for i=1:N
    bx(i)=bxc;
end
% Glen's coefficient in adjoint viscosity
nn=n;
% first point
i = 1;
ip = i + 1;
% 1st eq
A11(i,i)=-u(i)/dx;
A11(i,ip)=u(i)/dx;
A12(i,i)=-fic(i)/dx+bx(i)*rhoig;
A12(i,ip)=fic(i)/dx;
F1(i)=Fh(i);
% 2nd eq
A21(i,i)=H(i)/dx;
A21(i,ip)=-H(i)/dx;
etaHp=(eta(i+1)*H(i+1)+eta(i)*H(i))/(2*nn);
A22(i,i)=-2*etaHp/h2-beta(i)*m;
A22(i,ip)=etaHp/h2;
F2(i)=-Fu(i);
%    A(k,k)=1;
F2(i)=0;
% inner points
for i=2:N1
    ip = i + 1;
    im = i - 1;
    % 1st eq
    A11(i,im)=gam1(i)/dx;
    A11(i,i)=-u(i)/dx-2*gam1(i)/dx;
    A11(i,ip)=u(i)/dx+gam1(i)/dx;
    
    A12(i,im)=gam2(i)/dx;
    A12(i,i)=-fic(i)/dx+bx(i)*rhoig-2*gam2(i)/dx;
    A12(i,ip)=fic(i)/dx+gam2(i)/dx;
    F1(i)=Fh(i);
    % 2nd eq
    A21(i,im)=gam(i)/dx;
    A21(i,i)=H(i)/dx-2*gam(i)/dx;
    A21(i,ip)=-H(i)/dx+gam(i)/dx;
    etaHm=(eta(i-1)*H(i-1)+eta(i)*H(i))/(2*nn);
    etaHp=(eta(i+1)*H(i+1)+eta(i)*H(i))/(2*nn);
    A22(i,i)=-(etaHm+etaHp)/h2-beta(i)*m;
    A22(i,im)=etaHm/h2;
    A22(i,ip)=etaHp/h2;

    F2(i)=-Fu(i);
end
% last point
i = N;
im = i - 1;
% 1st eq
A11(i,i)=-u(i)/dx-2*gam1(i)/dx;
A12(i,im)=gam2(i)/dx;
A12(i,i)=-fic(i)/dx+bx(i)*rhoig-2*gam2(i)/dx;
F1(i)=Fh(i);

% 2nd eq
A21(i,im)=gam(i)/dx;
A21(i,i)=H(i)/dx-2*gam(i)/dx;
etaHm=(eta(i-1)*H(i-1)+eta(i)*H(i))/(2*n);
A22(i,i)=-2*etaHm/h2-beta(i)*m;
A22(i,im)=etaHm/h2;
F2(i)=-Fu(i);