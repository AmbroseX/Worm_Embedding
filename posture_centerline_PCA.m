%中心化
PX=[Xcenterline(:,:,1),Xcenterline(:,:,2)];
mm=mean(PX,1);   %对每一列进行平均
PX=PX-mm;

eignum=6;
%计算特征值，其中牵涉但数学变换
Cov=PX'*PX; %Cov n*n 因为n小一点 便于计算 
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

saveas(gcf, fullfile(savefolder,strcat(wormName,'_','posture','_','EigenValuePCA.jpg')));

% %显示特征worm
figure('Name',strcat(filepath,'_',wormName,'_','EigenWorm'),'NumberTitle','off');
for i=1:6
    si=num2str(i);
    x=V(1:100,i);
    y=V(101:end,i);
    subplot(2,3,i)
    plot(x,y)
    title(['EigenWormPCA' si]);
end
saveas(gcf, fullfile(savefolder,strcat(wormName,'_','posture','_','EigenWormPCA.jpg')));
