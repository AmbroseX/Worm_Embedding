function [y,Speed_downsampled_To_Energy]=Curve_DownSample_toPCA(curve,speed,v,kmode,downr)
% downr is ratio of folding
%输入原始的angledata，和PCA降维矩阵，进行PCA降维 投影

rtheta_s = curve; %100列维度
mean_X = mean(rtheta_s, 1); %整列做平均
T_sample = size(rtheta_s, 1);
X_shifted = rtheta_s - repmat(mean_X, T_sample, 1); %将rtheta_s减去均值，中心化

T = ceil(T_sample / downr); %downsampled data points
L = size(rtheta_s, 2); %number of segments 
X_downsampled = zeros(T, L);

Energy_speed=speed(:,2).*speed(:,2)+speed(:,3).*speed(:,3);

Speed_downsampled_To_Energy = zeros(T-1,1);
Speed_downsampled_To_Energy=decimate(Energy_speed,downr);

for i = 1:L
    X_downsampled(:, i) = decimate(X_shifted(:, i), downr); %decimate向下采样时会过滤掉再进行采样，防止混叠
end

%其元素的顺序颠倒。如果 A 是一个列向量，则 fliplr(A) 只返回 A，对于多维数组，fliplr
%作用于由第一个和第二个维度构成的平面。特征值由大到小的排列的向量
y = (X_downsampled * v);  %投影到eigenvectors张成的空间
y = y(:, 1:kmode); %the coefficients of the first K_modes

end