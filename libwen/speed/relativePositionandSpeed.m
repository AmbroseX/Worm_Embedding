function [relativePosion,relativeSpeed]=relativePositionandSpeed(lineposition,stageposition,times)
%算Centerline相对平台的位移
%A is wormdata
[m,n,~]=size(lineposition);

micronperunit = 0.05;  %每个单位位移代表的长度 um
micronperpxl = 2.43;  %每个像素代表的长度 um

lineposition=lineposition*micronperpxl;   %Centerline相对相机的时间-位移 %convert the pixel into micrometer

stageposition=stageposition*micronperunit;% wormdata.StagePosion

relativePosion=zeros(m,n,2);
relativeSpeed=zeros(m-1,n,2);
%calculate the actual relativePosition of the worm at each frame, the reference coordinate
%system here is as that of the moving stage

for i = 1:m
    relativePosion(i,:,1)=lineposition(i,:,1)-stageposition(i,1);
    relativePosion(i,:,2)=lineposition(i,:,2)-stageposition(i,2);
    if i>1 && i<m
        relativeSpeed(i-1,:,1)=(relativePosion(i,:,1)-relativePosion(i-1,:,1))/(times(i)-times(i-1));%calcute the speed,the time metric-s, the speed metric-μm/s
        relativeSpeed(i-1,:,2)=(relativePosion(i,:,2)-relativePosion(i-1,:,2))/(times(i)-times(i-1));
    end
end

end