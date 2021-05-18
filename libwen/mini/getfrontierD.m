function [flag,stotal]=getfrontierD(D,ratio)
num=length(D);
ss=0;
stotal=zeros(1,num);
for i=1:num
    ss=ss+D(i);
    if ss/sum(D)>ratio
        flag=i;
        break
    end
end

for i=1:num
    if i>1
        stotal(i)=stotal(i-1)+D(i);
    else
        stotal(i)=D(i);
    end   
end
stotal=stotal/sum(D);

end