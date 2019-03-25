clear
close all
%%
addpath('SSA')
%% Load final H and u from init file
load('DATA/SSAinit_N400.mat')
saveFlag = 0;

%%
pertubation = 0.01;
dC = pertubation.*C;

%% Setup restart
N_restart = 20000;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt = 1;
u_ref=u;

%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end

%%
plot(u_ref-u_mat(:,end))
%%
if saveFlag 
    save(['DATA/SSAPerturb_N', num2str(N),'_C' , num2str(pertubation*100,'%03.f') ,'.mat'], 'x','u','H','u_ref');
end
