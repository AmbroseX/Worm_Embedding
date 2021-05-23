%进行PCA降维，然后先进行了Embedding_K参数寻找，再通过Kstar寻找最佳mstar
%
savefolder=fullfile(workpath,'prodata','EmbeddingPhaseSpace');

km=str2num(wormName(2));
kmStart=[26,13;45,11;10,12;61,12];
%answer=inputdlg({'Kstar','mstar'},'PhaseSpace最佳嵌入参数');

Kstar=kmStart(km,1);
mStar =kmStart(km,2);  % to set best embedding dimension


downR=1;  %采样打几折 ratio of folding
lenr=1;   %画图打几折
 
%%% Our data 得到PCA投影矩阵
angle=angle_data;
ts=16; %sampling rate
rtheta_s = angle(:,2:end); %去掉第一列 得100列维度
mean_X = mean(rtheta_s, 1); %整列做平均
T_sample = size(rtheta_s, 1);
X_shifted = rtheta_s - repmat(mean_X, T_sample, 1); %将rtheta_s减去均值，中心化

K_mode = 6; % 取PCA前面 K_mode维度

T = ceil(T_sample /downR); %downsampled data points
L = size(rtheta_s, 2); %number of segments 
X_downsampled = zeros(T, L);

for i = 1:L
    X_downsampled(:, i) = decimate(X_shifted(:, i), downR); %decimate向下采样，每downR 个点取一个点
end
[eigenvectors,eigenvalues] = eig(X_downsampled' * X_downsampled);   %存储所有值特征向量是列向量，但是是从小达到排列的
eigenvectors = fliplr(eigenvectors); %得到从大到小排列的特征值矩阵，如果 A 是一个矩阵，反向排列
%V = (X_downsampled * eigenvectors);  %投影到eigenvectors张成的空间



[ybase,energy]=DownSample_toPCA(angle_data,speed,eigenvectors,K_mode,downR);  %PCA降维后的数据并取前面K_mode项

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


%% Get modes&trajectory by SVD 
%Using SVD to embedding
Ybase(isnan(Ybase)) = 0; %让Y中所有为nan的元素为0
[U, S, Vbase] = svd(Ybase, 'econ');  % V is the base



%导入计算相空间轨迹的其他angle_data数据
[wormfile,pathname]=uigetfile('.mat','选择要计算的文件',fullfile('workpath','data'));
testfilepath=strsplit(pathname,'\');  %
testfilepath=testfilepath{6};  %测试虫子种类

data=load(wormfile); % 1*12 cell ,33600*5 double
angle_X=data.wormdata.angle_data;
testName=data.wormdata.wormname;  %测试虫子编号
Xtime=data.wormdata.TimeElapsed;
[Xcenterline,Xspeed]=relativePositionandSpeed(data.wormdata.Centerline,data.wormdata.StagePosition,Xtime);
XFrame=data.wormdata.Framenum;
XcenterSpeed=centeroidSpeed(Xspeed,20); %质心的相对速度 um/s

[ytest,energy]=DownSample_toPCA(angle_X,Xspeed,eigenvectors,K_mode,downR);  %PCA降维后的数据并取前面K_mode项
%y是数据投影后的结果


%% Get Delay matrix by Wen
num_Frame = size(ytest, 1);   %y是前面K_modes个modes的矩阵系数
Ytest = zeros(num_Frame - K_delay + 1, K_mode_dimension * K_delay); %初始化test embedding后的矩阵大小

%embedding the data
j = 1;
for k = 1:K_delay
    Ytest(:, j:j + K_mode_dimension - 1) = ytest(k:k + num_Frame - K_delay, :);
    j = j + K_mode_dimension;
end
energyD=energy(K_delay:end);
ybartest=Ytest;  %embedding后的数据

FBtest=ForwardBackwardFrames(Xcenterline,Xtime,2);  %基线虫前进还是后退
FBtest=FBtest(K_delay:end-1);   %去除前面K_mode_dimension 个
XFrame=XFrame(K_delay:end-1);  %去除前面K_mode_dimension 个
SpeedFBtest=FB_Speed(Xcenterline,XcenterSpeed);
SpeedFBtest=SpeedFBtest(K_delay:end);

%% Get modes&trajectory by SVD
%Using SVD to embedding
Ytest(isnan(Ytest)) = 0; %让Y中所有为nan的元素为0

Xtestsvd = (Ytest * Vbase(:, 1:embedded_dimension))';


figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','X45'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
hold on
plot2d_color(Xtestsvd (4,:),Xtestsvd (5,:),SpeedFBtest)
xlabel('X4')
ylabel('X5')
hold off

saveas(gcf, fullfile(savefolder,'Angle_SVD',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5_line2','.jpg')))

% %plot
% xstatlen=4000;
% xendlen=10000;
% figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','X45'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
% hold on
% plot2d_color(Xtestsvd (4,xstatlen:xendlen),Xtestsvd (5,xstatlen:xendlen),SpeedFBtest(xstatlen:xendlen-1))
% xlabel('X4')
% ylabel('X5')
% xlim([-60,60])
% ylim([-60,60])
% hold off
% saveas(gcf, fullfile(savefolder,'Angle_SVD',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5_time_s',num2str(xstatlen),'_e',num2str(xendlen),'.jpg')))
% 
% 
% %plot time
% timelengh=3000;
% figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','X45_Time'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
% hold on
% plot2d_color_time(Xtestsvd (4,1:timelengh),Xtestsvd (5,1:timelengh),XFrame(1:timelengh-1))
% xlabel('X4')
% ylabel('X5')
% hold off
% 
% saveas(gcf, fullfile(savefolder,'Angle_SVD',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5_timeS_',num2str(timelengh),'.jpg')))

%plot X2
figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','X26'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
hold on
plot2d_color(Xtestsvd (2,:),Xtestsvd (6,:),SpeedFBtest)
xlabel('X2')
ylabel('X6')
hold off
saveas(gcf, fullfile(savefolder,'Angle_SVD',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X2_X6_line2','.jpg')))




%先画单条轨迹随着时间变化的图像
fignum=nchoosek(K_mode,2);  %图片数
if K_mode>mStar
   fignum=nchoosek(mStar,2);  %图片数
   K_mode=mStar;
end
len = floor(size(Xtestsvd , 2) / lenr)-1;
figcol=3; %三列
figrows=ceil(fignum/figcol); %行数
flag=1;
defig=figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','PhaseSpace'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
suptitle('2维投影X相空间轨迹(X降低大小后)');
hold on
for i=1:K_mode
    for j=1:K_mode
        if i<j
           if flag>fignum
              break
           end
           subplot(5,3,flag)
           plot(Xtestsvd (i, end - len:end), Xtestsvd (j, end - len:end));
           xlabel(strcat('X_',num2str(i)))
           ylabel(strcat('X_',num2str(j)))
           title(strcat('Xtestsvd -EigenWorm ',num2str(i)))
           flag=flag+1;
        end
    end
end
hold off

saveas(gcf, fullfile(savefolder,'Angle_SVD',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_PhaseSpace_line2','.jpg')))

close all

% 
% %先画单条轨迹随着时间变化的图像
% fignum=nchoosek(K_mode,2);  %图片数
% figcol=3; %三列
% figrows=ceil(fignum/figcol); %行数
% flag=1;
% nodefig=figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','PhaseSpace'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
% suptitle('2维投影X相空间轨迹');
% hold on
% for i=1:K_mode
%     for j=1:K_mode
%         if i<j
%            subplot(5,3,flag)
%            plot(Xtestsvd (i, :), Xtestsvd (j, :));
%            xlabel(strcat('X_',num2str(i)))
%            ylabel(strcat('X_',num2str(j)))
%            title(strcat('Xtestsvd -EigenWorm ',num2str(i)))
%            flag=flag+1;
%         end
%     end
% end
% hold off


% 
% len = floor(size(Xtestsvd , 2) / lenr)-1;
% 
% %先画单条轨迹随着时间变化的图像
% fignum=nchoosek(K_mode,2);  %图片数
% if K_mode>mStar
%    fignum=nchoosek(mStar,2);  %图片数
%    K_mode=mStar;
% end


% %绘制前进后退
% figcol=3; %三列
% figrows=ceil(fignum/figcol); %行数
% flag=1;
% defig=figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','PhaseSpace','_FB'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
% suptitle('2维投影X相空间轨迹(X降低大小后)');
% hold on
% for i=1:K_mode
%     for j=1:K_mode
%         if i<j
%            if flag>fignum
%               break
%            end
%            subplot(5,3,flag)
%            plot2d_FB_color(Xtestsvd (i, end - len:end), Xtestsvd (j, end - len:end),FBtest);
%            xlabel(strcat('X_',num2str(i)))
%            ylabel(strcat('X_',num2str(j)))
%            title(strcat('X-EigenWorm ',num2str(i)))
%            flag=flag+1;
%         end
%     end
% end
% hold off

%plot framenum dot

% plotframesFB
% 
% saveas(gcf, fullfile(savefolder,'FB_3',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Frame_FB','.jpg')))




% figure
% hold on
% num=length(FBtest);
% cmap =[0,0,1;0.3,1,0.6;1,0,0]; % colormap, with N colors  [1,0,0]红色,[0,0,1]蓝色,[0.3,1,0.6]绿色
% for i = 1:num
%     if FBtest(i)>0
%         plot(Xtestsvd (4,i),Xtestsvd (5,i),'.','color',cmap(3,:))
%     elseif FBtest(i)==0
%         plot(Xtestsvd (4,i),Xtestsvd (5,i),'.','color',cmap(2,:))
%     else
%         plot(Xtestsvd (4,i),Xtestsvd (5,i),'.','color',cmap(1,:))
%     end
% end
% xlabel('X4')
% ylabel('X5')
% colormap(cmap)
% h=colorbar;
% set(get(h,'label'),'string','Backward-->Stand or Turn-->Forward');
% hold off

% saveas(gcf, fullfile(savefolder,'FB_3',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5_dot','.jpg')))
% 
% %plot line
% figure
% hold on
% num=length(FBtest);
% cmap =[0,0,1;0.3,1,0.6;1,0,0]; % colormap, with N colors  [1,0,0]红色,[0,0,1]蓝色,[0.3,1,0.6]绿色
% plot2d_FB_color(Xtestsvd (4,:),Xtestsvd (5,:),FBtest)
% xlabel('X4')
% ylabel('X5')
% colormap(cmap)
% h=colorbar;
% set(get(h,'label'),'string','Backward-->Stand or Turn-->Forward');
% hold off
% 
% saveas(gcf, fullfile(savefolder,'FB_3',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5_line','.jpg')))
% 
% 
% 
% 
% 
% 
% 
% 
% figure
% hold on
% plot2d_FB_color(Xtestsvd (4,:),Xtestsvd (5,:),FBtest)
% xlabel('X4')
% ylabel('X5')
% hold off
% 
% saveas(gcf, fullfile(savefolder,'FB',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X4_X5','.jpg')))
% 
% figure
% hold on
% plot2d_FB_color(Xtestsvd (2,:),Xtestsvd (6,:),FBtest)
% xlabel('X4')
% ylabel('X5')
% hold off
% 
% saveas(gcf, fullfile(savefolder,'FB',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X2_X6','.jpg')))
% 
% figure
% hold on
% plot2d_FB_color(Xtestsvd (5,:),Xtestsvd (6,:),FBtest)
% xlabel('X5')
% ylabel('X6')
% hold off
% 
% saveas(gcf, fullfile(savefolder,'FB',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_X5_X6','.jpg')))
% 
% 
% figure
% hold on
% title('SVD phase plot')
% plot(Xtestsvd (1, end - len:end), Xtestsvd (2, end - len:end), Xtestsvd (3, end - len:end));
% xlabel('X_t_1')
% ylabel('X_t_2')
% zlabel('X_f_2')
% hold off
% saveas(gcf,fullfile(savefolder,'SVD_phase.jpg'))



% 
% figure
% hold on
% plot2d_color(Xtestsvd (2,:),Xtestsvd (6,:),energyD)
% hold off


