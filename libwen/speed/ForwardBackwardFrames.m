function y=ForwardBackwardFrames(centerline,time,fr)
%给定centerline的数据,time,和用几帧的数据来计算
% fr 前后的帧数
colnum=10;  %前面多少个点来进行计算
[m,~,~]=size(centerline);
y=zeros(m,1);%存储前进还是后退
for i=1:m
    if i>fr && i<=m-fr
        y_head=ForwardBackwardFrame(centerline(i-fr:i+fr,1:colnum,:),time(i-fr:i+fr),fr); %head 
        y_tail=ForwardBackwardFrame(centerline(i-fr:i+fr,end-colnum:end,:),time(i-fr:i+fr),fr);  %tail
        if y_head>0 && y_tail>0  %头尾都前进就前进
            y(i)=1;  
        elseif y_head<0 && y_tail <0 %头尾都后退才算后退
            y(i)=-1;
        else          % 算作不动或者做turn
            y(i)=0;
        end
    else
        if rand(1)>0.5  %前面几帧随机选取
            y(i)=1;
        else 
            y(i)=-1;
        end      
    end
end

  
end



