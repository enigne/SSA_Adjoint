% close all
clear

%% Set range
C_range = [10:1:1030];
b_range = [10:1:1030];
u_range = [10:1033];
h_range = u_range;
saveFlag = 1;
N = 1600;

%% Load data
load('DATA/SSA_Macayeal_Adjoint_N1600_u.mat')
wuc = wght_mat(u_range, C_range);
wuc_dcouple = vWght_mat(u_range, C_range);

load('DATA/SSA_Macayeal_Adjoint_N1600_H.mat')
whc = wght_mat(h_range, C_range);

W_mat = {wuc,  wuc_dcouple};
nameList = {'\widetilde{\mathbf{\Sigma}}_{uC}', '\widehat{\mathbf{\Sigma}}_{uC}'};


%% SVD
ni = length(W_mat);
xC = length(C_range);
xu = length(u_range);

S_mat = zeros(xC, ni);
U_mat = zeros(xu, ni);
V_mat = zeros(xC, ni);

for i=1: ni
    [U, S, V] = svd(W_mat{i});
    U_mat(:,i) = U(:, end);
	S_mat(:,i) = diag(S);
    V_mat(:,i) = V(:, end);
end

%%
if saveFlag
    save(['DATA/SSA_Macayeal_Adjoint_N', num2str(N), '_S_SVD.mat'], 'S_mat', ...
        'U_mat', 'V_mat', 'nameList');
end
%%
figure
semilogy(S_mat)
legend(nameList)
ylim([1e-12,1e-3])
xlim([0,1024])
% figure
% subplot(1, 2, 1)
% % semilogy(abs(U_mat),'o')
% plot((U_mat),'o')
% title('U')
% legend(nameList, 'Location', 'best')
% subplot(1, 2, 2)
% % semilogy(abs(V_mat),'o')
% plot((V_mat),'o')
% title('V')
% legend(nameList, 'Location', 'best')
% 
% 
%%