clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N1600.mat')

%% Hack for m=1
C = beta;
m = 1;
%% Compute analytical solutions
saveFlag = 1;
iGL = GLpos;
slicePos = 900;
% [u_ana,H_ana,beta_ana,x_ana,phi_ana,psi_ana,wght_ana,bwght_ana,dudC]=...
%     Hum(N,L,H(2),H(iGL+1),iGL,slicePos,C,rhog,as,1.00, m);

[x_ana, H_ana, u_ana, psi_ana, phi_ana, wght_ana, bwght_ana] = ...
    analyticalSSA(H(2), H(iGL+1), m, C(2:end-1), as, rhog, x(2:end-1), iGL, slicePos);
%% plot
figure
subplot(2,1,1)
plot(x, u)
hold on
plot(x_ana, u_ana)
subplot(2,1,2)
plot(x, H)
hold on
plot(x_ana, H_ana)
figure
subplot(3,1,1)
plot(x_ana, psi_ana)
subplot(3,1,2)
plot(x_ana, phi_ana)
subplot(3,1,3)
plot(x_ana,wght_ana) 
%% save
if saveFlag
    save(['DATA/SSAAnalytical_N', num2str(N), '.mat'], ...
         'x_ana', 'u_ana', 'H_ana', 'phi_ana', 'psi_ana', 'wght_ana', 'bwght_ana', 'slicePos');
end