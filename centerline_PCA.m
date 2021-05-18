

% To get the centerline EigenWorm
numfram=length(centerline);
X=zeros(numfram,200); %存储 BoundaryA X, BoundaryB X,BoundaryA Y,BoundaryB Y
for i=1:numfram
    X(i,:)=[centerline(i,:,1),centerline(i,:,2)]; %M*400
end
[m,n]=size(X);
eignum=200;

%中心化
mm=mean(X,1);   %对每一列进行平均
Xmean=repmat(mm,m,1); %M*200
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
saveas(gcf, fullfile(savefolder,strcat(wormName,'_','centerline','_','EigenValuePCA.jpg')));

save(fullfile(savefolder,strcat(wormName,'_','centerline_PCA.mat')),'V'); % %存储矩阵

% %显示特征line
figure('Name',strcat(filepath,'_',wormName,'_','EigenWorm'),'NumberTitle','off');
for i=1:6
    si=num2str(i);
    x=V(1:100,i);
    y=V(101:end,i);
    subplot(2,3,i)
    plot(x,y)
    title(['EigenWormPCA' si]);
end
saveas(gcf, fullfile(savefolder,strcat(wormName,'_','centerline','_','EigenWormPCA.jpg')));