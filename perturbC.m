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
lWin = 900e3;
rWin = 1000e3;
C = C.*abs(u).^(m-1);
dC = pertubation.*C.*((x>=lWin)&(x<=rWin));

%% Pertubation in time
N_restart = 300;
dt_pert = 0.05;
sInd = [0: (N_restart-1)];
seasonType = 0;

if seasonType == 0
    % sine [-1,1]
    season = sin(sInd*dt_pert*2*pi);
elseif seasonType == 1
    % cosine [-1,1]
    season = cos(sInd*dt_pert*2*pi);
else
    % square [0, 1]
    season = (sin(sInd*dt_pert*2*pi) > 0)*2-1;
end

%% Setup restart
H_mat = zeros(length(H), N_restart);
u_mat = zeros(length(u), N_restart);
gpos_vec = zeros(1, N_restart);
u_ref= u;
H_ref = H;
%% Solve SSA GL problem
for i = 1: N_restart
    [gpos, H, u, beta]=FlowlineSSA(H, b, x, dx, Nx, A, C+dC*season(i), 1, n, rhoi, ...
        rhow, g, as, dt_pert, dt_pert, u);
    H_mat(:, i) = H;
    u_mat(:, i) = u;
    gpos_vec(i) = gpos;
end

%% Plot
tMesh = linspace(0, N_restart*dt_pert, N_restart);
[Xm, Tm] = meshgrid(x, tMesh);
u_diff = u_mat - u_ref;
H_diff = H_mat - H_ref;
figure
subplot(2, 1, 1)
contour(Xm, Tm, u_diff', 100)
colorbar
subplot(2, 1, 2)
contour(Xm, Tm, H_diff', 100)
colorbar

%%
figure
plot(x, u_mat(:,end) -u_ref)

%% Cut data
% cutInd = [1:ceil(10/dt)+1, N_restart];
% u_mat = u_mat(:,cutInd);
% H_mat = H_mat(:,cutInd);
% 
%% Save
if saveFlag 
    save(['DATA/SSAPerturb_N', num2str(N),'_C' , num2str(pertubation*100,'%03.f'),...
        '_x', num2str(lWin/1000) ,'.mat'], 'x','u_mat','H_mat','u_ref','H_ref',...
        'dt_pert','N_restart','pertubation', 'dC', 'season');
end
