%saveas(gcf, fullfile(savefolder,'ICA_Curve',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_PhaseSpace','.jpg')))

%导入计算相空间轨迹的其他angle_data数据
[wormfile,pathname]=uigetfile('.mat','选择要计算的文件',fullfile('workpath','data'));
testfilepath=strsplit(pathname,'\');  %
testfilepath=testfilepath{6};  %测试虫子种类

data=load(wormfile); % 1*12 cell ,33600*5 double
curve_X=data.wormdata.curv_data;
testName=data.wormdata.wormname;  %测试虫子编号
Xtime=data.wormdata.TimeElapsed;
[Xcenterline,Xspeed]=relativePositionandSpeed(data.wormdata.Centerline,data.wormdata.StagePosition,Xtime);
XcenterSpeed=centeroidSpeed(Xspeed,20); %质心的相对速度 um/s
XFrame=data.wormdata.Framenum;
[ytest,energy]=Curve_DownSample_toPCA(curve_X,Xspeed,eigenvectors,K_mode,downR);  %PCA降维后的数据并取前面K_mode项
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

Ztestica = deMix * Ytest' + (deMix * mu') * ones(1, size(Ytest, 1));
save(fullfile(workpath,'prodata','curve_ICA',strcat(filepath,'_',wormName,'_','curve_ICAbaseZ.mat')),'Ztestica'); % %存储矩阵


len = floor(size(Ztestica , 2) / lenr)-1;
%先画单条轨迹随着时间变化的图像
fignum=nchoosek(K_mode,2);  %图片数
figcol=3; %三列
figrows=ceil(fignum/figcol); %行数
flag=1;
nodefig=figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','PhaseSpace'),'NumberTitle','off');  %基线虫种类_编号_测试线虫种类_编号
suptitle('2维投影X相空间轨迹');
hold on
for i=1:K_mode
    for j=1:K_mode
        if i<j
           subplot(figrows,3,flag)
           plot(Ztestica (i, end-len:end), Ztestica (j, end-len:end));
           xlabel(strcat('Z_',num2str(i)))
           ylabel(strcat('Z_',num2str(j)))
           title(strcat('Ztestica -EigenWorm ',num2str(i)))
           flag=flag+1;
        end
    end
end
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_PhaseSpace','.jpg')))


figure
hold on
plot2d_color(Ztestica (1,:),Ztestica (2,:),SpeedFBtest)
xlabel('Z1')
ylabel('Z2')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z1_Z2','.jpg')))

figure
hold on
plot2d_color(Ztestica (1,:),Ztestica (3,:),SpeedFBtest)
xlabel('Z1')
ylabel('Z3')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z1_Z3','.jpg')))

figure
hold on
plot2d_color(Ztestica (3,:),Ztestica (4,:),SpeedFBtest)
xlabel('Z3')
ylabel('Z4')
hold off

figure
hold on
plot2d_color(Ztestica (3,:),Ztestica (5,:),SpeedFBtest)
xlabel('Z3')
ylabel('Z5')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z3_Z4','.jpg')))


figure
hold on
plot2d_color(Ztestica (1,:),Ztestica (5,:),SpeedFBtest)
xlabel('Z1')
ylabel('Z5')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z1_Z5','.jpg')))

figure
hold on
plot2d_color(Ztestica (1,:),Ztestica (6,:),SpeedFBtest)
xlabel('Z1')
ylabel('Z6')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z1_Z6','.jpg')))

figure
hold on
plot2d_color(Ztestica (2,:),Ztestica (3,:),SpeedFBtest)
xlabel('Z2')
ylabel('Z3')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z2_Z3','.jpg')))

figure
hold on
plot2d_color(Ztestica (2,:),Ztestica (5,:),SpeedFBtest)
xlabel('Z2')
ylabel('Z5')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z2_Z5','.jpg')))

figure
hold on
plot2d_color(Ztestica (2,:),Ztestica (6,:),SpeedFBtest)
xlabel('Z2')
ylabel('Z6')
hold off

saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_Z2_Z6','.jpg')))


figure
hold on
plot2d_color(Ztestica (4,:),Ztestica (5,:),SpeedFBtest)
xlabel('Z4')
ylabel('Z5')
hold off



figure
hold on
[Kymo_curvdata, Kymo_theta] = Get_Kymograph(mix, Kymoratio, eigenvectors, K_delay, K_mode, embedded_dimension);  %得到模式图
hold off
saveas(gcf, fullfile(savefolder,'curve_ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_mode','.jpg')))