answer=inputdlg('请输入最佳Kstar','The best mstar');
mStar = 6;  % to set best embedding dimension

%center ybar by removing mean
ybar = bsxfun(@minus, ybar, mean(ybar));
Y = ybar;

%% Get Delay matrix by Wen
Kymoratio_Jiahao = 100;
K_delay = 12;
K_mode_dimension = 5;
embedded_dimension = mStar;
num_Frame = size(y, 1);   %y是前面K_modes个modes的矩阵系数
Y = zeros(num_Frame - K_delay + 1, K_mode_dimension * K_delay); %初始化embedding后的矩阵大小

%embedding the data
j = 1;
for k = 1:K_delay
    Y(:, j:j + K_mode_dimension - 1) = y(k:k + num_Frame - K_delay, :);
    j = j + K_mode_dimension;
end
ybar=Y;


%% Get modes&trajectory by SVD or ICA
%using covairance matrix
%[V,D] = eig(Y'*Y);
%V = fliplr(V);
%c_embedded = Y*V;

Y(isnan(Y)) = 0; %让Y中所有为nan的元素为0
[U, S, V] = svd(Y, 'econ');
XX = (Y * V(:, 1:embedded_dimension))';
len = floor(size(XX, 2) / 2);


figure
hold on
title('SVD phase plot')
plot3(XX(1, end - len:end), XX(2, end - len:end), XX(4, end - len:end));
xlabel('X_t_1')
ylabel('X_t_2')
zlabel('X_f_2')
hold off
saveas(gcf,fullfile(savefolder,'SVD_phase.jpg'))