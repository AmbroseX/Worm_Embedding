close all
tic
%% -------------load data-------------------------------
clear all
clc
% N2 silence testDLP


%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embedding');
[wormfile,pathname]=uigetfile('.mat','选择要计算的文件',fullfile('workpath','data'));

filepath=strsplit(pathname,'\');
filepath=filepath{6};

%For the 2080Ti
% workpath=fullfile('/','home','wenlab','xrk','Worm_Embed');

addpath(genpath(fullfile(workpath,'libwen')));
addpath(genpath(fullfile(workpath,'data')));
addpath(genpath(pathname));

disp('Staring load data...')

load(wormfile) % 1*12 cell ,33600*5 double

wormName = wormdata.wormname;   %to create folder to keep .jpg
savefolder=fullfile(workpath,'prodata',filepath,wormName);

if exist(savefolder)==0
    disp('dir is not exist');
    mkdir(savefolder);
    disp('make dir success');
else
    disp('dir is exist');
end
% Our data: preparation
angle_data=wormdata.angle_data;
curve_data=wormdata.curv_data;
time=wormdata.TimeElapsed;
Frame=wormdata.Framenum;
[centerline,speed]=relativePositionandSpeed(wormdata.Centerline,wormdata.StagePosition,time);  %算出来绝对速度
centerSpeed=centeroidSpeed(speed,25);

disp('loadding success ...')
%% -------------------load success---------------------------------

%% get the best K and m

%posture_PCA
% angle_PCA;


%Embedding_K;
%saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_','Tpred(K).jpg')))   %for angle
%saveas(gcf, fullfile(workpath,'prodata','curve',strcat(filepath,'_',wormName,'_','Tpred(K).jpg')))   % for curve


%using the code to calculate best m
%Embedding_m;
% saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_','Tpred(m).jpg')))
% saveas(gcf, fullfile(workpath,'prodata','curve',strcat(filepath,'_',wormName,'_','Tpred(m).jpg')))   % for curve


%Embedding_PhaseSpace
%saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_XPahseSpace_do','_',num2str(downR),'_de_',num2str(lenr),'.jpg')))  %de是减小了大小的意思
%saveas(gcf, fullfile(savefolder,'FB',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_XPahseSpace_FB.jpg')))  %de是减小了大小的意思

% Embedding_PhaseSpace_ICA
% saveas(gcf, fullfile(savefolder,'ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_ZPahseSpace.jpg')))
% [Kymo_curvdata, Kymo_theta] = Get_Kymograph(mix, Kymoratio, eigenvectors, K_delay, K_mode, embedded_dimension);  %得到模式图
% saveas(gcf, fullfile(savefolder,'ICA',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_ZMode.jpg')))


%Curve_Embedding_PhaseSpace

%Curve_Embedding_PhaseSpace_ICA