%进行PCA降维，然后先进行了Embedding_K参数寻找，再通过Kstar寻找最佳mstar
%
savefolder=fullfile(workpath,'prodata','EmbeddingPhaseSpace');

km=str2num(wormName(2));
kmStart=[12,13;45,11;10,12;61,12];
%answer=inputdlg({'Kstar','mstar'},'PhaseSpace最佳嵌入参数');

Kstar=kmStart(km,1);
mStar =kmStart(km,2);  % to set best embedding dimension


downR=1;  %采样打几折 ratio of folding
lenr=2;   %画图打几折
K_mode =5; % 取PCA前面 K_mode维度


%%% Our data 得到PCA投影矩阵
ts=16; %sampling rate
rtheta_s = curve_data; %100列维度
mean_X = mean(rtheta_s, 1); %整列做平均
T_sample = size(rtheta_s, 1);
X_shifted = rtheta_s - repmat(mean_X, T_sample, 1); %将rtheta_s减去均值，中心化



T = ceil(T_sample /downR); %downsampled data points
L = size(rtheta_s, 2); %number of segments 
X_downsampled = zeros(T, L);

for i = 1:L
    X_downsampled(:, i) = decimate(X_shifted(:, i), downR); %decimate向下采样，每downR 个点取一个点
end
[eigenvectors,eigenvalues] = eig(X_downsampled' * X_downsampled);   %存储所有值特征向量是列向量，但是是从小达到排列的
eigenvectors = fliplr(eigenvectors); %得到从大到小排列的特征值矩阵，如果 A 是一个矩阵，反向排列
%V = (X_downsampled * eigenvectors);  %投影到eigenvectors张成的空间



[ybase,energy]=Curve_DownSample_toPCA(curve_data,speed,eigenvectors,K_mode,downR);  %PCA降维后的数据并取前面K_mode项

%% Get Delay matrix by Wen
K_delay = Kstar;
K_mode_dimension = K_mode;
embedded_dimension = mStar;

num_Frame = size(ybase, 1);   %ybase是基的前面K_modes个PCA modes的矩阵系数
Ybase = zeros(num_Frame - K_delay + 1, K_mode_dimension * K_delay); %初始化embedding后的矩阵大小

%embedding the data
j = 1;
for k = 1:K_delay
    Ybase(:, j:j + K_mode_dimension - 1) = ybase(k:k + num_Frame - K_delay, :);
    j = j + K_mode_dimension;
end
energyD=energy(K_delay:end);
ybarbase=Ybase;  %embedding后的数据

FBbase=ForwardBackwardFrames(centerline,time,2);  %基线虫前进还是后退
FBbase=FBbase(K_delay:end-1);   %去除前面K_mode_dimension 个
Framebase=Frame(K_delay:end-1);  %去除前面K_mode_dimension 个
SpeedFBbase=FB_Speed(centerline,centerSpeed);
SpeedFBbase=SpeedFBbase(K_delay:end);




%% ICA by our method:Centered and whitened
%mu = mean(X',2);
%T = sqrtm(inv(cov(X)));
%Zcw = T * bsxfun(@minus,X',mu);
Kymoratio=100;
mu = mean(Ybase);
[Xbase, whiteningMatrix, dewhiteningMatrix] = whitenRows(bsxfun(@minus, Ybase, mu)', embedded_dimension);
[~, Bbase] = fastICA(Xbase, embedded_dimension, 'negentropy');
deMix = Bbase' * whiteningMatrix;  %白化矩阵
mix = dewhiteningMatrix * Bbase;
Zbaseica = deMix * Ybase' + (deMix * mu') * ones(1, size(Ybase, 1));

ICAbase.mix=mix;
ICAbase.mu=mu;
save(fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_','curve_ICAbase.mat')),'ICAbase'); % %存储矩阵
save(fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_','curve_ICAbaseZ.mat')),'Zbaseica'); % %存储矩阵

len = floor(size(Zbaseica , 2) / lenr)-1;
%画base相空间轨迹随着时间变化的图像
fignum=nchoosek(K_mode,2);  %图片数
figcol=3; %三列
figrows=ceil(fignum/figcol); %行数
flag=1;
nodefig=figure('Name',strcat('Base',filepath,'_',wormName,'_','_','ICA_PhaseSpace'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
suptitle('2维投影X相空间轨迹');
hold on
for i=1:K_mode
    for j=1:K_mode
        if i<j
           subplot(figrows,3,flag)
           plot(Zbaseica (i, end-len:end), Zbaseica (j, end-len:end));
           xlabel(strcat('Z_',num2str(i)))
           ylabel(strcat('Z_',num2str(j)))
           title(strcat('Zbaseica -EigenWorm ',num2str(i)))
           flag=flag+1;
        end
    end
end
hold off

saveas(gcf, fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_PhaseSpace','.jpg')))


figure
hold on
plot2d_color(Zbaseica (1,:),Zbaseica(5,:),SpeedFBbase)
xlabel('Z1')
ylabel('Z5')
hold off

saveas(gcf, fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_Z1_Z5','.jpg')))



figure
hold on
plot2d_color(Zbaseica (3,:),Zbaseica(4,:),SpeedFBbase)
xlabel('Z3')
ylabel('Z4')
hold off

saveas(gcf, fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_Z3_Z4','.jpg')))

figure
hold on
plot2d_color(Zbaseica (2,:),Zbaseica(3,:),SpeedFBbase)
xlabel('Z2')
ylabel('Z5')
hold off

saveas(gcf, fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_Z2_Z5','.jpg')))


figure
hold on
plot2d_color(Zbaseica (4,:),Zbaseica(5,:),SpeedFBbase)
xlabel('Z4')
ylabel('Z5')
hold off



figure
hold on
plot2d_color(Zbaseica (3,:),Zbaseica(5,:),SpeedFBbase)
xlabel('Z3')
ylabel('Z5')
hold off

figure
hold on
[Kymo_curvdata, Kymo_theta] = Get_Kymograph(mix, Kymoratio, eigenvectors, K_delay, K_mode, embedded_dimension);  %得到模式图
hold off
saveas(gcf, fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_mode','.jpg')))