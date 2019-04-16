function [x, H, u, psi, v, wght, bwght] = analyticalSSA(H0, HGL, m, C, a, rhog, x, iGL, ist)
%% prepare for \int_x^x* Cx^{m-1} dx
GLmask = (x < x(iGL+1));
integrant = C.*x.^m .* GLmask;
Hfact = trapz(x, integrant) - cumtrapz(x, integrant);

%% Compute for H
H =  (HGL^(m+2)  + (m+2)*a^m/rhog * Hfact).^(1/(m+2));

%% u
u = a*x./H;

%% psi
xst = x(ist+1);
Hst = H(ist+1);
mask = (x > xst);

fic = xst/(rhog*Hst^(m+3));
axstm = (a*xst)^m;
psi = -1/Hst+C.*fic.*((a.*x).^m-axstm);

% integrant = C.*x.^(m-1) .* (1-mask);
% psifact = trapz(x, integrant) - cumtrapz(x, integrant);
% psi = -1/Hst - m*a^(m+1).* fic.*psifact;

psi(mask) = 0;

%% v
v = -a*fic*H.^m;
v(mask) = 0;
%% C weights
wght =- v .* u.^m;

%% b weights
bwght = (m+1)*C*a*xst.*(a*x).^m./(rhog*Hst^(m+3).*H);

