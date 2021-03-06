function [A11, A12, A21, A22, F1, F2, ux, eta, beta]=constrauctAdjSSAMatrices...
    (N, n, ist, sigma, u, H, bxc, Afact, rhoig, dx, uObs, HObs, glInd, epsilon, Cbeta, m, MacyealFlag, transientFlag)
if nargin < 18
    transientFlag = 1;
    if nargin < 17
        MacyealFlag = 0;
    end
end
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
Afact=2/Afact^(1/n);
mm1=m-1;
ux=zeros(N,1);
eta=ux;
beta=ux;
fic=ux;
bx=ux;
F = ux;
%
p=(1-n)/n;
N1=N-1;
hh=2*dx;
% u derivative
for i=2:N1
    ux(i)=(u(i+1)-u(i-1))/hh;
end
ux(1)=ux(2);
ux(N)=ux(N1);

% gaussian for the observation
gaussfact=1/sqrt(2*pi*sigma^2);
xi = zeros(N, 1);
for i=1:N
    xi(i)=(i-ist)*dx;
    F(i)=gaussfact*exp(-xi(i)^2/(2*sigma^2));
end
Fint = dx*sum(F);
F = F./Fint;
% check which observation function
if uObs
    Fu = F;
else
    Fu = 0*F;
end
if HObs
    Fh = F;
else
    Fh = 0*F;
end

% viscosity, fi coeff in first eq, bottom derivative
eta = Afact*abs(ux).^p;
fic = eta .* ux -rhoig*H;
beta = Cbeta .* abs(u).^mm1;
bx = bxc* ones(N,1);

% add artifical viscosity
% epsilon
% D2(N,dx)
visco = zeros(N,1);
if transientFlag
    viscosInd = 80;
else
    viscosInd = 0;
end
for i = 1:N
    if abs(i - glInd) < viscosInd
        visco(i) = epsilon;
    end
end
ArtiV = dx * visco .* D2(N,dx)  ;

% construct matrix
if MacyealFlag
    Heta = 1./1 .* H .* eta;
else
    Heta = 1./ n .* H .* eta;
end
dHeta = Dcd(N,dx)*Heta;
A11 = u .* Dup(N, dx, u) + ArtiV;
A12 = fic .* Dup(N, dx, fic) + spvardiag(bx)*rhoig;

A21 = - H .* Dup(N, dx, H);
A22 = Heta .* (Dp(N, dx) * Dm(N,dx)) + dHeta.* Dup(N,dx,dHeta) - m * spvardiag(beta) + ArtiV;
F1=Fh;
F2=-Fu;


%% Boundary conditions
% psi(N) = 0
A11(N,:) = 0;
A11(N, N) = A11(N-1,N-1);
A12(N,:) = 0;
F1(N) = 0;

% phi(0)'= psi(0)' = 0
A11(1, 2) = -A11(1, 1);
A12(1, 2) = -A12(1, 1);
A21(1, 2) = -A21(1, 1);
A22(1, 2) = -A22(1, 1);
F2(1) = 0;
F1(1) = 0;

% phi(N) = 0
A21(N, :) = 0;
A22(N, :) = 0;
A22(N, N) = A22(N-1, N-1);
F2(N) = 0;
