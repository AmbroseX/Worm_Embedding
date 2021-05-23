function y=FB_postureFrames(centerline,fr)
%给定centerline的数据,time,和用几帧的数据来计算
% fr 前后的帧数
fr=1;
colnum=10;  %前面多少个点来进行计算
[m,~,~]=size(centerline);
y=zeros(m,1);%存储前进还是后退
for i=1:m
    if i>fr && i<=m-fr
    y(i)=FB_postureFrame(centerline(i-fr:i+fr,:,:),colnum);
    else
        if rand(1)>0.5  %前面几帧随机选取
            y(i)=1;
        else 
            y(i)=-1;
        end
    end
end

frames=wormdata.Framenum;
figure
hold on
for i=1:length(y)
    if y(i)>0
        plot(frames(i),0.1*y(i),'r.');  
    else
        plot(frames(i),0.1*y(i),'b.');
    end
end
xlim([10000,14000])
ylim([-2,2])
hold off
    
    
end





end



