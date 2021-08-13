x=[0 1 1 0];
y=[0 0 1 1];  %定义四个点 [0 0] [1 0] [1 1] [0 1]
H_F=fill(x,y,[0 0.1 0.2 0.6]);  %定义四个点的C值

row_cmap = 100;  %定义色图矩阵的行数
color_map1=zeros(row_cmap,3);  %定义色图矩阵
color_r = 0:1/(row_cmap-1):1; 
color_g = 0:1/(row_cmap-1):1;
color_b = 0:1/(row_cmap-1):1;
color_map1(:,1) = color_r; 
color_map1(:,3) = color_g;
colormap(color_map1);
colorbar;

t=[linspace(0,2*pi) nan];
x=sin(t);y=cos(2*t);z=sqrt(t);%所要绘制的曲线方程
patch(x,x,y,z,'edgecolor','flat','facecolor','none')
colormap(jet)
xlabel('x')
ylabel('y')
zlabel('z')
view(3);
grid on;
colorbar;


x = -10:0.1:10;
y = x.^2;
patch(x,y,y,'EdgeColor','flat','LineWidth',1,'MarkerFaceColor','flat','FaceColor','none')
grid on; colorbar






n = 100; 
x = linspace(-10,10,n); y = x.^2; 
p = plot(x,y,'r', 'LineWidth',5); 
% modified jet-colormap 
cd = [uint8(jet(n)*255) uint8(ones(n,1))].'; 
drawnow 
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd) 





%% 以下是一种可能的方法：使用从所需颜色图中获取的不同颜色明确绘制线的每个线段。
x = -10:10; % x data. Assumed to be increasing 
num = length(x); % number of colors. Assumed to be greater than size of x 
y = x.^2; % y data
z=x;
c = exp(x)-1;  % speed;
%c(end)=[];
% cc储存，第一列数值，第二列其原本序号，第3列颜色序号
cc =zeros(num,3);
cc(:,1)=c;
cc(:,2)=1:num;
cc=sortrows(cc,1);  %按照第1列排序，第二列是其原本序号序列
cc(:,3)=1:num;   %第三列得到颜色排序序列
cc=sortrows(cc,2);  %再按照第2列排序回来

cmap = jet(num ); % colormap, with N colors
linewidth = 0.5; % desired linewidth 
xi = x(1)+linspace(0,1,num +1)*x(end); % interpolated x values 
yi = interp1(x,y,xi); % interpolated y values  所有的插值方法都要求x是单调的，并且xi不能够超过x的范围。
hold on 
for n = 1:num-1 
    plot(x([n n+1]), y([n n+1]), 'color', cmap(cc(n,3),:), 'linewidth', linewidth); 
end
colormap(cmap)
h=colorbar;
set(get(h,'title'),'string','mm/s')
caxis([-1,1])
hold off


plot3(x,y,z)
plot3d_color(x,y,z,c);

hold on
fig1=plot2d_color(x,y,c);
xlabel('x')
hold off

x=-10:10;
y=0:20;
num=length(x);
cmap=jet(30);
cmap=[0,0,1;1,0,0];
linewidth = 0.5; % desired linewidth 
figure
hold on 
for n = 1:num-1
    fig=plot(x([n n+1]), y([n n+1]), 'color', cmap(mod(n,2)+1,:), 'linewidth', linewidth); 
end
colormap(cmap)
h=colorbar;
set(get(h,'label'),'string','Backward->Forward');



function popupmenu1_Callback(hObject, eventdata, handles)
val=get(hObject,'value');%获取数值,从上到下依次1到4
str=get(hObject,'String');%获取字符串，这里是菜单所有的字符串，相当于存到了字符串数组里
switch val
    case 1
        set(handles.edit1,'String','');%输入到可编辑文本里
    case 2
        set(handles.edit1,'String',str{2});
    case 3
        set(handles.edit1,'String',str{3});
    case 4
        set(handles.edit1,'String',str{4});
end


figure
hold on
x=-20:5;
y=x.*x;
v=x;
v(end)=[];
plot2d_color(x,y,v)
hold off



