%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embedding');
addpath(genpath(fullfile(workpath,'libwen')));
savefolder=fullfile('G:','Data','WenLab','Worm_Embedding','prodata','PhaseSpace');

filepath={'N2','silence'};
worm={'w1','w2','w3','w4'};
sumN2=3;
sumSilence=5;
%load N2
for i=1:sumN2
    N2_PCA_data{i}=load(fullfile(workpath,'prodata',filepath{1},worm{i},strcat(filepath{1},'_',worm{i},'_','angle_PCA.mat')));
end
for i=1:sumN2
    PCA_matrix{i}=N2_PCA_data{i}.V;
end

clear N2_PCA_data

mean_PCA_matrix=PCA_matrix{1};
for i=2:sumN2
    mean_PCA_matrix=mean_PCA_matrix+PCA_matrix{i};
end
mean_PCA_matrix=mean_PCA_matrix/sumN2;

clear PCA_matrix

wormName='N2';
% %显示N2 的平均 worm
figure('Name','N2:平均EigenWorm','NumberTitle','off');
for i=1:6
    si=num2str(i);
    x=linspace(0,1,100);   %0-之间的100个点
    y=mean_PCA_matrix(:,i);
    subplot(2,3,i)
    plot(x,y)
    title(['EigenWormPCA ' si]);
end
saveas(gcf, fullfile(savefolder,strcat('N2','_','mean_angle','_','EigenWormPCA.jpg')));


%load N2 原始数据
workfolder=fullfile(workpath,'data','N2');
filename=dir(fullfile(workfolder,'*.mat'));
for i=1:sumN2
   dataN2{i}=load(fullfile(workfolder,filename(i).name));
end

for i=1:sumN2
    angle_dataN2{i}=dataN2{i}.wormdata.angle_data(:,2:end); %去掉第一列 得100列维度
    time_elapsedN2{i}=dataN2{i}.wormdata.TimeElapsed;
    speed_dataN2{i}=dataN2{i}.wormdata.speed;  %time,Vx,Vy
end


%load silence 原始数据
workfolder=fullfile(workpath,'data','silence');
filename=dir(fullfile(workfolder,'*.mat'));
for i=1:sumSilence
    silence_PCA_data{i}=load(fullfile(workfolder,filename(i).name));
end
for i=1:sumSilence
    angle_dataSilence{i}=silence_PCA_data{i}.wormdata.angle_data(:,2:end); %去掉第一列 得100列维度
    time_elapsedSilence{i}=silence_PCA_data{i}.wormdata.TimeElapsed;
    speed_dataSilence{i}=silence_PCA_data{i}.wormdata.speed;  %time,Vx,Vy
end








