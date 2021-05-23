function fig=plot2d_FB_color(x,y,fb)
%将前进和后退点映射为颜色图，前进红色，后退蓝色
% x data. Num个数据
% y data  Num个数据
% v 根据这个数值画颜色图，有 Num-1 个
num = length(x); % number of colors. Assumed to be greater than size of x


cmap =[0,0,1;0.3,1,0.6;1,0,0]; % colormap, with N colors  [1,0,0]红色,[0,0,1]蓝色,[0.3,1,0.6]绿色
linewidth = 0.5; % desired linewidth 
hold on 
for n = 1:num-1
    if fb(n)>0  %前进
        r=3;
    elseif fb(n)==0  %原地不动或者转弯
        r=2;
    else   %后退
        r=1; 
    end
    fig=plot(x([n n+1]), y([n n+1]), 'color', cmap(r,:), 'linewidth', linewidth); 
end
colormap(cmap)
h=colorbar;
set(get(h,'label'),'string','Backward-->Stand or Turn-->Forward');
hold off

end
