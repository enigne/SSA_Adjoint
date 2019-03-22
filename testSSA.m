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
gpos_vec = zeros(1, N_restart);

%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end
%%


% For adjssa.m you need the input
rhoig = rhoi*g;
u = u_mat(:,1);
H = H_mat(:,1);
sigma=0.5e4; 
ist = 200;
n=3;

x = x(2:end);
u = u(2:end);
% H on stagger grid
H = (H(1:end-1)+H(2:end)) * 0.5;
beta = beta(2:end);

%%

[psi,fi,wght,bwght,A,f,fic,eta,beta,fix,psix,hx,ux,B,G,K]=adjssa(Nx-1,n,ist,sigma,u,H,mean(bxc),A,rhoig,dx);
