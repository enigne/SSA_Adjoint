clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')

%% Setup restart
N_restart = 100;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
%% Solve SSA GL problem
for i = 1: N_restart
    [~, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
end
%%
% %%
% x0=0; 
% xN=1.8e6; 
% h=dx; 
% rhoig=9.8*900; 
% AA=1.4648e-17; 
% bxc=constant b_x; 
% sigma=0.5e4; 
% xst=x(ist); 
% n=3; 
% iGL=123;
% 
% [psi,fi,wght,bwght,A,f,fic,eta,beta,fix,psix,hx,ux,B,G,K]=adjssa(N,n,ist,xst,sigma,u,H,betain,bxc,AA,rhoig,h,x0,xN);
