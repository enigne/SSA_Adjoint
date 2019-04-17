function [x, H, u, psi, v, wght, bwght] = analyticalSSA(H0, HGL, m, C, a, rhog, x, iGL, ist)
%% prepare for \int_x^x* Cx^{m-1} dx
xGL = x(iGL+1);
GLmask = (x < xGL);
integrant = C.*x.^m .* GLmask;
Hfact = trapz(x, integrant) - cumtrapz(x, integrant);

%% Compute for H
H =  (HGL^(m+2)  + (m+2)*a^m/rhog * Hfact).^(1/(m+2));

%% u
u = a*x./H;

%% psi
xst = x(ist+1);
Hst = H(ist+1);
mask = (x < xst);

fic = xst/(rhog*Hst^(m+3));
axstm = (a*xGL)^m;
psi =  C.*fic.*(axstm - (a.*x).^m);
% 
% integrant = C.*x.^(m-1).*GLmask;
% psifact =  trapz(x, integrant) - cumtrapz(x, integrant);
% psi = m*a^(m+1).* fic.*psifact;

psi(mask) = -1/Hst + psi(ist+1);

%% v
v = a*fic*H.^m;
v(mask) = 0;
%% C weights
wght =- v .* u.^m;

%% b weights
bwght = (m+1)*C*a*xst.*(a*x).^m./(rhog*Hst^(m+3).*H);

