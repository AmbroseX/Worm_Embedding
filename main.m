close all
tic
%% -------------load data-------------------------------
clear,clc
% N2 silence testDLP


%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embedding');
[wormfile,pathname]=uigetfile('.mat','选择要计算的文件',fullfile('workpath','data'));

filepath=strsplit(pathname,'\');
filepath=filepath{6};

%For the 2080Ti
% workpath=fullfile('/','home','wenlab','xrk','Worm_Embed');

addpath(genpath(fullfile(workpath,'libwen')));
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
speed=wormdata.speed;
centerline=wormdata.Centerline;
disp('loadding success ...')
%% -------------------load success---------------------------------

%% get the best K and m

%posture_PCA
%Embedding_K;
% saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_','Tpred(K).jpg')))
angle_PCA;