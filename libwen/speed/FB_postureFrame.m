function FB=FB_postureFrame(center,col)
%i是第i时刻
%利用前后各1帧或者两帧来计算此时是前进还是后退
%fr是利用前后一帧或者两帧来计算
[m,~,~]=size(center);

position=zeros(4,2);  %存储 t-2时刻的Head，Tail位置，  Neck 的t-2和t+2时刻位置
position(1,:)=center(1,1,:);  % head0
position(2,:)=center(1,end,:);  % tail0
position(3,:)=center(1,col,:);   % neck0
position(4,:)=center(m,col,:);  % neck1

d_neck1_head0=norm(position(4,:)-position(1,:));
d_neck1_tail0=norm(position(4,:)-position(2,:));
d_neck0_head0=norm(position(3,:)-position(1,:));
d_neck0_tail0=norm(position(3,:)-position(2,:));

if d_neck1_head0<d_neck0_head0 && d_neck1_tail0>d_neck0_tail0
    FB=1;
else
    FB=-1;
end
end