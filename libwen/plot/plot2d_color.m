function fig=plot2d_color(x,y,v)
% x data. Num个数据
% y data  Num个数据
% v 根据这个数值画颜色图，有 Num-1 个
num = length(x); % number of colors. Assumed to be greater than size of x
nv1=sum(v>0);
nv2=length(v)-nv1;

%disp(length(x))
%disp(length(y))
% cc储存，第一列速度数值，第二列其原本序号，第3列颜色序号
cc =zeros(num-1,3);
cc(:,1)=v;% speed;归一化
cc(:,2)=1:num-1;
cc=sortrows(cc,1);  %按照第1列升序排序，第二列是其原本序号序列
cc(:,3)=1:num-1;   %第三列得到颜色排序序列
cc=sortrows(cc,2);  %再按照第2列排序回来，用第三列颜色编号来绘图

cmap = jet(num-1); % colormap, with N colors
linewidth = 1; % desired linewidth 
%xi = x(1)+linspace(0,1,num +1)*x(end); % interpolated x values 
%yi = interp1(x,y,xi); % interpolated y values  所有的插值方法都要求x是单调的，并且xi不能够超过x的范围。

hold on 
for n = 1:num-1
    fig=plot(x([n n+1]), y([n n+1]), 'color', cmap(cc(n,3),:), 'linewidth', linewidth); 
end
colormap(cmap)
h=colorbar;
set(get(h,'title'),'string','mm/s')
caxis([1e-3*min(v),1e-3*max(v)])
hold off

end
