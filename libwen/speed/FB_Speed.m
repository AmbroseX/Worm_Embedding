function FBspeed=FB_Speed(centerline,centerspeed)
colnum=10;  %前面多少个点来进行计算
[mc,~,~]=size(centerline);

FBspeed=zeros(mc-1,1);

xvec=[centerline(:,1,1)-centerline(:,colnum,1),centerline(:,1,2)-centerline(:,colnum,2)];  %第colnum点指向Head的向量

for i=1:mc-1
   FBspeed(i)=dot(xvec(i,:),centerspeed(i,:))/norm(dot(xvec(i,:),centerspeed(i,:)))*norm(centerspeed(i,:));  % dot(x,v)/abs(x*v)*abs(v)
end

end