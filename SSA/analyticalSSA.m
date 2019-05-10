function [x, H, u, psi_u, v_u, wght_u, bwght_u, psi_h, v_h, wght_h, bwght_h] = analyticalSSA(H0, HGL, m, C, a, rhog, x, iGL, ist)
%% prepare for \int_x^x* Cx^{m-1} dx
xGL = x(iGL+1);
GLmask = (x < xGL);
integrant = C.*x.^m .* GLmask;
Hfact = trapz(x, integrant) - cumtrapz(x, integrant);
dx = abs(x(2) - x(1));

%% Compute for H
H =  (HGL^(m+2)  + (m+2)*a^m/rhog * Hfact).^(1/(m+2));

%% u
u = a*x./H;

%% Masks
xst = x(ist+1);
Hst = H(ist+1);
maskST = (x < xst);
maskGL = (x > xGL);
mask = maskST | maskGL;

% 
integrant = C.*x.^(m-1).*GLmask;
integrant(1) = 0;
psifact =  trapz(x, integrant) - cumtrapz(x, integrant);

%% v
Cv = a/rhog.*xst./Hst.^(m+3);
v_u = Cv.*H.^m;
v_u(mask) = 0;

%% psi
psi_u = Cv .* m .* a^(m-1) .* psifact;
psi_u(maskST) = -1/Hst + psi_u(ist+1);
psi_u(maskGL) = 0;

%% C weights
wght_u =- v_u .* u.^m;

%% b weights
Hx = -C.*a.^m/rhog.*x.^m./H.^(m+1);
vx = m*Cv.*H.^(m-1).*Hx;
vx(maskST) = 0;
vx(ist+1) = u(ist)/rhog/Hst^2/dx;
bwght_u = rhog*(vx.*H + Hx.*v_u);
%% for h-response
%% v
Cv_h = -1./rhog/(Hst^(m+1));
v_h = Cv_h.*H.^m;
v_h(mask) = 0;

%% psi
psi_h =  Cv_h .* m .* a^(m-1) .* psifact;
psi_h(maskST) = psi_h(ist+1);
psi_h(maskGL) = 0;

%% C weights
wght_h =- v_h .* u.^m;

%% b weights
Hx = -C.*a.^m/rhog.*x.^m./H.^(m+1);
vx_h = m*Cv_h.*H.^(m-1).*Hx;
vx_h(maskST) = 0;
% vx_h(ist+1) = -1./rhog/Hst/dx;
bwght_h = rhog*(vx_h.*H + Hx.*v_h);
% bwght_h(ist+1) = bwght_h(ist+1) + 1/dx;
