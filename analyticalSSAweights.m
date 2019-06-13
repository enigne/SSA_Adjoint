%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute SSA perturbations du and dH for given db and dC
%
% Author: Cheng Gong
% Date: 2019-06-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, du, dH] = analyticalSSAweights(HGL, m, C, a, rhog, x, iGL, A, n, db, dC)
%% prepare for \int_x^x* Cx^{m-1} dx
xGL = x(iGL+1);
GLmask = (x < xGL);
intWeights = C.*x.^m .* GLmask;
Hfact = trapz(x, intWeights) - cumtrapz(x, intWeights);
dx = abs(x(2) - x(1));

%% Compute for H
H =  (HGL^(m+2)  + (m+2)*a^m/rhog * Hfact).^(1/(m+2));

%% u
u = a*x./H;

%% prepare for integrate weights
intWeights = C.*(a*x).^(m).*((m+1).*db./H + dC./(C+1e-16));
intWeights = intWeights .* GLmask;
weights =  (trapz(x, intWeights) - cumtrapz(x, intWeights))./(rhog*H.^(m+1));

%% eta is needed for dh
Hx = -C.*a.^m/rhog.*x.^m./H.^(m+1);
% correction at ist for psi
eta = 2 * A^(-1/n).*(a.*(H-x.*Hx)./H.^2).^((1-n)/n);
% D1 operator
Dx = Dup(length(x), dx, u);

%% compute du and dH
du = u./H.*db -u./H.*weights;
dH = eta./(n*rhog.*H).*(Dx*(u.*db)/dx/dx).*x + weights;
