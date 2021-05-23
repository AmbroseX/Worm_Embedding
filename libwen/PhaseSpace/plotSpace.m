function figcolor=plotSpace(trajectory,PCA_matrix,speed,wormname,w,savefolder,savename)
    %savefolder is 存储的路径 savename是angle_data，还是Embedding_data
    speed=sqrt(speed(:,2).*speed(:,2)+speed(:,3).*speed(:,3)); %转换为速度的大小
    maxspeed=max(speed);
    speed=speed/maxspeed; %归一化
    %startPointColor = [0,0,1];  %蓝色
    %endPointColor = [1,0,0];  %红色
    %sz=5; %size of point
    
    project=trajectory*PCA_matrix;
    
%     figure
%     hold on
%     plot3d_color(project(:,1),project(:,2),project(:,3),speed)
%     hold off
    
    %绘图
    figure
    hold on
    plot2d_color(project(:,1),project(:,2),speed)
    xlabel('PCA-1')
    ylabel('PCA-2')
    title(wormname)
    hold off
    saveas(gcf, fullfile(savefolder,strcat('PhaseSpace','_',savename,'_',wormname,'_',w,'_','12.jpg')))
    
    
    figure
    hold on
    plot2d_color(project(:,1),project(:,3),speed);
    xlabel('PCA-1')
    ylabel('PCA-3')
    title(wormname)
    hold off
    saveas(gcf, fullfile(savefolder,strcat('PhaseSpace','_',savename,'_',wormname,'_',w,'_','13.jpg')))
    
    figure
    hold on
    plot2d_color(project(:,2),project(:,3),speed);
    xlabel('PCA-2')
    ylabel('PCA-3')
    title(wormname)
    hold off
    saveas(gcf, fullfile(savefolder,strcat('PhaseSpace','_',savename,'_',wormname,'_','23.jpg')))
    
    project=project(1:end-1,:); %去掉最后一行画3d
    
    figure
    hold on
    patch(project(:,1),project(:,2),project(:,3),speed,'edgecolor','flat','facecolor','none');
    view(3);
    xlabel('PCA-1')
    ylabel('PCA-2')
    zlabel('PCA-3')
    title(wormname)
    grid on;
    colormap(jet);
    colorbar;
    hold off
    
    
    saveas(gcf, fullfile(savefolder,strcat('PhaseSpace','_',savename,'_',wormname,'_',w,'_','123.jpg')))
%     figcolor{1}=fig1;
%     figcolor{2}=fig2;
%     figcolor{3}=fig3;
%     figcolor{4}=fig4;
end