% To get the angle EigenWorm
angle=angle_data;
rtheta_s = angle(:,2:end); %去掉第一列 得100列维度
numfram=length(rtheta_s);

X=rtheta_s; %存储 M*100


[m,n]=size(X);
eignum=100;

%中心化
mm=mean(X,1);   %对每一列进行平均
Xmean=repmat(mm,m,1); %M*100
X=X-Xmean;


%计算特征值，其中牵涉但数学变换
Cov=X'*X; %Cov n*n 因为n小一点 便于计算 
[V,D]=eigs(Cov,eignum); %最大的eignum个 V是特征向量400*400，D是特征值400*1
D=diag(D);
[flag,s]=getfrontierD(D,0.95);

%绘制前95%的eigenvalue
figure('Name','EigenValue','NumberTitle','off');
hold on
plot(s(1:flag+1))
title(strcat(filepath,'-',wormName,'-','worm'))
xlabel('Eigennum')
ylabel('Sum of EigenValue')
hold off
saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_','angle','_','EigenValuePCA.jpg')));

save(fullfile(savefolder,strcat(filepath,'_',wormName,'_','angle_PCA.mat')),'V'); % %存储矩阵

% %显示特征worm
figure('Name',strcat(filepath,'_',wormName,'_','EigenWorm'),'NumberTitle','off');
for i=1:6
    si=num2str(i);
    x=linspace(0,1,100);   %0-之间的100个点
    y=V(:,i);
    subplot(2,3,i)
    plot(x,y)
    title(['EigenWormPCA' si]);
end
saveas(gcf, fullfile(savefolder,strcat(filepath,'_',wormName,'_','angle','_','EigenWormPCA.jpg')));