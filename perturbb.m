clear
close all
%%
addpath('SSA')
N = 1600;
%% Load final H and u from init file
load(['DATA/SSAinit_N', num2str(N), '.mat'])
saveFlag = 1;

%% Perturbation on b 0.01meter
pertubation = 0.01;
lWin = 900e3;
rWin = 1000e3;
% use linear sliding law with variable C(x)
C = C.*abs(u).^(m-1);
% db = pertubation.*((x>=lWin)&(x<=rWin)).*(1-abs(x-9.5e5)/5e4);
db = pertubation.*((x>=lWin)&(x<=rWin));

%% Setup restart
N_restart = 30000;
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
dt_pert = 0.5;
u_ref= u;
%% unperturb H
% !!! This is very important for transient simulation since we are
% acctually looking at h, the thickness should be adjusted to have h as it
% was from the initialization, otherwise sensitivity analysis can not
% predict the behavior. However, for steady state, this does not matter.
H = H - db;
H_ref = H;

%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b+db, x, dx, Nx, A, C, 1, n, rhoi, ...
        rhow, g, as, dt_pert, dt_pert, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end

%% Plot
figure
subplot(2,1,1)
plot(x, u_mat(:,2) - u_ref)
subplot(2,1,2)
plot(x, H_mat(:,2) - H_ref)

%% Cut data
cutInd = [1:ceil(10/dt_pert)+1, N_restart];
u_mat = u_mat(:,cutInd);
H_mat = H_mat(:,cutInd);

%% Save
if saveFlag 
    save(['DATA/SSAPerturb_N', num2str(N),'_b' , num2str(pertubation*100,'%03.f'),...
        '_x', num2str(lWin/1000) ,'.mat'], 'x','u_mat','H_mat','u_ref','H_ref',...
        'dt_pert','N_restart','pertubation', 'db');
end
