function centerSpeed=centeroidSpeed(Speed,binbar)
if nargin <2,binbar=1;end
%质心的时间-位移 X1
mx=Speed(:,:,1);
my=Speed(:,:,2);
centerSpeed(:,1)=smooth(mean(mx,2),binbar); % 行平均
centerSpeed(:,2)=smooth(mean(my,2),binbar); %列平均
%做smooth平滑处理
end