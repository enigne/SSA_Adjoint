clear
close all
%%
addpath('SSA')
N = 1600;
%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 1;

%%
pertubation = 0.01;
lWin = 900e3;
rWin = 1000e3;
C = C.*abs(u).^(m-1);
dC = pertubation.*C.*((x>=lWin)&(x<=rWin));

%% Setup restart
N_restart = 300;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt_pert = 0.05;
u_ref= u;
H_ref = H;
%% Solve SSA GL problem
for i = 1: N_restart
    % Seasonal variation
    season = sin(i*dt_pert*2*pi);
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC*season, 1, n, rhoi, ...
        rhow, g, as, dt_pert, dt_pert, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end

%% Save
if saveFlag 
    save(['DATA/SSASeasonal_N', num2str(N),'_C' , num2str(pertubation*100,'%03.f'),...
        '_x', num2str(lWin/1000) ,'.mat'], 'x','u_mat','H_mat','u_ref','H_ref',...
        'dt_pert','N_restart','pertubation', 'dC');
end
