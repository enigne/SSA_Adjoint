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
saveFlag = 0;
iGL = GLpos;
slicePos = 700;

[x_ana, H_ana, u_ana, psi_ana, phi_ana, wght_ana, bwght_ana, psi_ana_h, phi_ana_h, wght_ana_h, bwght_ana_h] = ...
    analyticalSSA(H(end-1), H(iGL+1), m, C(2:end-1), as, rhog, x(2:end-1), iGL, slicePos, A, n, rhow, rhoi);

%% Plot solutions
figure
subplot(2,1,1)
plot(x, u)
hold on
plot(x_ana, u_ana)
subplot(2,1,2)
plot(x, H)
hold on
plot(x_ana, H_ana)
% figure
% subplot(4,1,1)
% plot(x_ana, psi_ana)
% subplot(4,1,2)
% plot(x_ana, phi_ana)
% subplot(4,1,3)
% plot(x_ana,psi_ana_h)
% subplot(4,1,4)
% plot(x_ana,phi_ana_h)
% figure
% subplot(4,1,1)
% plot(x_ana, wght_ana)
% subplot(4,1,2)
% plot(x_ana, bwght_ana)
% subplot(4,1,3)
% plot(x_ana,wght_ana_h)
% subplot(4,1,4)
% plot(x_ana,bwght_ana_h)
% 
% %% Pertubations
% pertubation = 0.01;
% lWin = 900e3;
% rWin = 1000e3;
% dC = pertubation.*C.*((x>=lWin)&(x<=rWin));
% db = pertubation.*((x>=lWin)&(x<=rWin));
% 
% %% Compute analytical solutions
% [~, dudC_ana, dHdC_ana] = analyticalSSAweights( H(iGL+1), m, C(2:end-1), ...
%         as, rhog, x(2:end-1), iGL, A, n, 0, dC(2:end-1));
%     
% [~, dudb_ana, dHdb_ana] = analyticalSSAweights( H(iGL+1), m, C(2:end-1), ...
%         as, rhog, x(2:end-1), iGL, A, n, db(2:end-1), 0);
% 
% %%
% figure
% subplot(2,1,1)
% plot(x_ana, dudC_ana)
% subplot(2,1,2)
% plot(x_ana, dHdC_ana)
% 
% figure
% subplot(2,1,1)
% plot(x_ana, dudb_ana)
% subplot(2,1,2)
% plot(x_ana, dHdb_ana)
% 
% 
% %% save
% if saveFlag
%     save(['DATA/SSAAnalytical_N', num2str(N), '.mat'], ...
%         'x_ana', 'u_ana', 'H_ana', 'phi_ana', 'psi_ana', 'wght_ana', 'bwght_ana', ...
%         'phi_ana_h', 'psi_ana_h', 'wght_ana_h', 'bwght_ana_h', 'slicePos', ...
%         'dudC_ana', 'dHdC_ana', 'dudb_ana', 'dHdb_ana');
% end