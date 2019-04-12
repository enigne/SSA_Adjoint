clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N1600.mat')

%% Hack for m=1
C = C.*abs(u).^(m-1);
m = 1;
%% Compute analytical solutions
saveFlag = 1;
iGL = GLpos;
slicePos = 500;
[u_ana,H_ana,beta_ana,x_ana,phi_ana,psi_ana,wght_ana,bwght_ana,dudC]=Hum(N,L,H(2),H(iGL),iGL,slicePos,C,rhog,as,1.001, m);

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
plot(x_ana,psi_ana)
figure
plot(x_ana,wght_ana) 
%% save
if saveFlag
    save(['DATA/SSAAnalytical_N', num2str(N), '.mat'], ...
         'x_ana', 'u_ana', 'H_ana', 'phi_ana', 'psi_ana', 'wght_ana', 'bwght_ana', 'slicePos');
end