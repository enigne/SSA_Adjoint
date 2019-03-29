clear
close all
%%
addpath('SSA')
N = 800;
%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 1;

%%
pertubation = 0.01;
lWin = 0e3;
rWin = 1600e3;
dC = pertubation.*C.*((x>=lWin)&(x<=rWin));

%% Setup restart
N_restart = 20000;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt = 1;
u_ref= u;
H_ref = H;
%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC, m, n, rhoi, rhow, g, as, dt, dt, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end

%%
plot(x, u_mat(:,end) -u_ref)
%%
if saveFlag 
    save(['DATA/SSAPerturb_N', num2str(N),'_C' , num2str(pertubation*100,'%03.f'), '_x', num2str(lWin/1000) ,'.mat'], 'x','u_mat','H_mat','u_ref','H_ref', 'dt','N_restart','pertubation', 'dC');
end
