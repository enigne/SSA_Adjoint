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
saveFlag = 1;
iGL = GLpos;
slicePos = [250, 500, 700, 900];
psi_a_mat = zeros(N, length(slicePos)*2);
phi_a_mat = zeros(N, length(slicePos)*2);
wght_a_mat = zeros(N, length(slicePos)*2);
bwght_a_mat = zeros(N, length(slicePos)*2);


%%
for i = 1:length(slicePos)
    [x_ana, H_ana, u_ana, psi_ana, phi_ana, wght_ana, bwght_ana, psi_ana_h, phi_ana_h, wght_ana_h, bwght_ana_h] = ...
        analyticalSSA(H(2), H(iGL+1), m, C(2:end-1), as, rhog, x(2:end-1), iGL, slicePos(i));

    %% scale to MPa unit
    psi_a_mat(:,i*2-1) = psi_ana;
    phi_a_mat(:,i*2-1) = 1e6*phi_ana;
    wght_a_mat(:,i*2-1) = 1e6*wght_ana;
    bwght_a_mat(:,i*2-1) = 1e6*bwght_ana;

    psi_a_mat(:,i*2) = psi_ana_h;
    phi_a_mat(:,i*2) = 1e6*phi_ana_h;
    wght_a_mat(:,i*2) = 1e6*wght_ana_h;
    bwght_a_mat(:,i*2) = 1e6*bwght_ana_h;
    
    %% plot
%     figure(1)
%     subplot(1, 2, 1)
%     plot(x_ana, psi_ana)
%     hold on
%     xlim([0,1.04e6])    
%     subplot(2, 1, 2)
%     plot(x_ana,psi_ana_h) 
%     hold on
%     xlim([0,1.04e6])    
%     
%     figure(2)
%     subplot(2, 1, 1)
%     plot(x_ana, phi_ana)
%     hold on
%     xlim([0,1.04e6])    
%     subplot(2, 1, 2)
%     plot(x_ana,phi_ana_h) 
%     hold on
%     xlim([0,1.04e6])    
% 
%     figure(3)
%     subplot(2, 1, 1)
%     plot(x_ana, wght_ana)
%     hold on
%     xlim([0,1.04e6])    
%     subplot(2, 1, 2)
%     plot(x_ana,wght_ana_h) 
%     hold on
%     xlim([0,1.04e6])    
% 
%     figure(4)
%     subplot(2, 1, 1)
%     plot(x_ana, bwght_ana)
%     hold on
%     xlim([0,1.04e6])    
%     subplot(2, 1, 2)
%     plot(x_ana,bwght_ana_h) 
%     hold on
%     xlim([0,1.04e6])    

end

%% save
if saveFlag
    save(['DATA/SSAAnalyticalMat_N', num2str(N), '.mat'], ...
         'x_ana', 'u_ana', 'H_ana', 'psi_a_mat', 'phi_a_mat', ...
         'wght_a_mat', 'bwght_a_mat', 'slicePos');
end
