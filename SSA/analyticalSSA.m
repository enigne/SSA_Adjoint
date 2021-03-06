function [x, H, u, psi_u, v_u, wght_u, bwght_u, psi_h, v_h, wght_h, bwght_h] = analyticalSSA(Hcf, HGL, m, C, a, rhog, x, iGL, ist, A, n, rhow, rhoi)
%% prepare for \int_x^x* Cx^{m-1} dx
xGL = x(iGL+1);
GLmask = (x < xGL);
integrant = C.*x.^m .* GLmask;
Hfact = trapz(x, integrant) - cumtrapz(x, integrant);
dx = abs(x(2) - x(1));
% dx = 1;

%% Compute for H
H =  (HGL^(m+2)  + (m+2)*a^m/rhog * Hfact).^(1/(m+2));

%% u
u = a*x./H;

%% Masks
xst = x(ist+1);
Hst = H(ist+1);
maskST = (x < xst);
maskGL = (x >= xGL);
mask = maskST | maskGL;

%% add floating and calving conditions
HGL = HGL *0.8;
H(maskGL) = HGL + (Hcf - HGL)/(x(end)-xGL)*(x(maskGL)-xGL);
% u(maskGL) = u(iGL+1) + A*(rhog/4*(1-rhoi/rhow))^n * (1/(n+1)*(x(end)-xGL)/(Hcf-HGL))*...
%         ((HGL+(Hcf - HGL)/(x(end)-xGL)*(x(maskGL)-xGL)).^(n+1) - HGL^(n+1));

% H(maskGL) = Hcf(maskGL);
u(maskGL) = A*(rhog/4*(1-rhoi/rhow))^n * cumtrapz(x(maskGL),H(maskGL).^n)+u(iGL+1);

%% 
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
% correction at ist for psi
eta_ist = 2 *A^(-1/n)*(a*(Hst-xst*Hx(ist+1))/Hst^2)^((1-n)/n);
psi_h(ist+1) = psi_h(ist+1) - eta_ist/n/rhog/Hst/dx;
%% C weights
wght_h =- v_h .* u.^m;

%% b weights
Hx = -C.*a.^m/rhog.*x.^m./H.^(m+1);
vx_h = m*Cv_h.*H.^(m-1).*Hx;
vx_h(maskST) = 0;
% vx_h(ist+1) = -1./rhog/Hst/dx;
bwght_h = rhog*(vx_h.*H + Hx.*v_h);
%
delta = eta_ist/n/rhog/Hst*u(ist+1);
bwght_h(ist) = bwght_h(ist) - delta/dx/dx;
bwght_h(ist+2) = bwght_h(ist+2) + delta/dx/dx;
