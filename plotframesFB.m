frames=wormdata.Framenum;
FB=ForwardBackwardFrames(centerline,wormdata.TimeElapsed,2);
%figure('Name',strcat(filepath,'_',wormName,'_',testfilepath,'_',testName,'_','Forward or Backward'),'NumberTitle','off'); 
hold on
for i=1:length(FB)
    if FB(i)>0
        plot(frames(i),0.1*FB(i),'r.');  
    elseif FB(i)==0
        plot(frames(i),0,'g.')
    else
        plot(frames(i),0.1*FB(i),'b.');
    end
end
xlabel('Frame')
ylabel('FB')
legend('Backward','Stand or Turn','Foward')
%xlim([8900,9200])
ylim([-2,2])
hold off


%三种状态占比
Fratio=sum(FB==1)/length(FB);
Bratio=sum(FB==-1)/length(FB);
Sratio=sum(FB==0)/length(FB);
