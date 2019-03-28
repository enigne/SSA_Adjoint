clear
% close all
%%
addpath('SSA')
%%
N = 800;
uObs = 1;
transientFlag = 0;
N_restart = 1;
pertubation = 0.01;
lWin = 0;

%%
if uObs
    obsName = '_u';
else
    obsName = '_H';
end

if transientFlag
    tName =  ['_T', num2str(N_restart)];
else
    tName = '';
end

pertName = ['_C' , num2str(pertubation*100,'%03.f')];
if lWin > 0
    lWinName = ['_x', num2str(lWin)];
else
    lWinName = '_x0';
end

% forward Pertubation
load(['DATA/SSAPerturb_N', num2str(N), pertName, lWinName, '.mat'])
% adjoint
load(['DATA/SSAAdjoint_N', num2str(N), tName, obsName, '.mat'])
% Initialization
load(['DATA/SSAinit_N', num2str(N), '.mat'])
%%
sensitivityC = dx*wght_mat(:,:,end)'*dC(2:end);

if uObs
    sol_ref = u_ref;
    sol_mat = u_mat;
else 
    sol_ref = H_ref;
    sol_mat = H_mat;
end

if transientFlag
    sol = sol_mat(:, 1);
else
    sol = sol_mat(:, end);
end
%%
figure
% hold on
plot(x, sol - sol_ref, 'LineWidth',2)
hold on
plot(x(ist), sensitivityC,'ro')
xlim([min(x), max(x)])
%%
% tempW = wght_mat(2:540,2:540);
% [U,S,V]=svd(tempW);
% figure(7);
% semilogy(diag(S),'bo')
% hold on
% ylim([1e-10,1e-4])