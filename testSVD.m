close all
clear

%% Set range
C_range = [10:1:1030];
b_range = [10:1:1030];
u_range = [10:1033];
h_range = u_range;
saveFlag = 0;
N = 1600;

%% Load data
load('DATA/SSAAdjoint_N1600_u.mat')
wub = bwght_mat(u_range, b_range);
wuc = wght_mat(u_range, C_range);

load('DATA/SSAAdjoint_N1600_H.mat')
whb = bwght_mat(h_range, b_range);
whc = wght_mat(h_range, C_range);

W_mat = {wuc, whc, wub, whb};
nameList = {'uC', 'hC', 'ub', 'hb'};
%%
A = [wub, wuc;
    whb, whc];

[U, S, V] = svd(A);

%%
figure
semilogy(diag(S));
figure
subplot(2,2,1)
semilogy(abs(U(:,end)),'o')
title('U')

subplot(2,2,2)
semilogy(abs(V(:,end)),'o')
title('V')

subplot(2,2,3)
plot(U(:,end),'o')
title('U')

subplot(2,2,4)
plot(V(:,end),'o')
title('V')

%%

% %% SVD
% ni = length(W_mat);
% xC = length(C_range);
% xu = length(u_range);
% 
% S_mat = zeros(xC, ni);
% U_mat = zeros(xu, ni);
% V_mat = zeros(xC, ni);
% 
% for i=1: ni
%     [U, S, V] = svd(W_mat{i});
%     U_mat(:,i) = U(:, end);
% 	S_mat(:,i) = diag(S);
%     V_mat(:,i) = V(:, end);
% end
% 
% %%
% if saveFlag
%     save(['DATA/SSAAdjoint_N', num2str(N), '_S_SVD.mat'], 'S_mat', ...
%         'U_mat', 'V_mat', 'nameList');
% end
% %%
% figure
% semilogy(S_mat)
% legend(nameList)
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
% %%