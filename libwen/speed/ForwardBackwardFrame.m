function y=ForwardBackwardFrame(center,t,fr)
%i是第i时刻
%利用前后各1帧或者两帧来计算此时是前进还是后退
%fr是利用前后一帧或者两帧来计算
[m,n,~]=size(center);
vecline=zeros(m-1,n,2);

for i=1:2*fr
    vecline(i,:,1)=(center(i+1,:,1)-center(i,:,1))/(t(i+1)-t(i));
    vecline(i,:,2)=(center(i+1,:,2)-center(i,:,2))/(t(i+1)-t(i));
end

%vx=vecline(:,:,1);vy=vecline(:,:,2);
xvec_to_head_0=[center(1,1,1)-center(1,end,1),center(1,1,2)-center(1,end,2)]; %开始帧的到头的向量
v=mean(vecline,2);  %这前面n个点的平均速度
v=reshape(v,[],2);
v=mean(v,1); %这几帧的平均速度

FB=dot(v,xvec_to_head_0);  %内积
if FB>0   %内积大于零是前进
    y=1;
else     %内积小于零是后退
    y=-1;
end
end
