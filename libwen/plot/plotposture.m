

%centerline=wormdata.Centerline;
x=centerline(:,:,1);
y=centerline(:,:,2);

figure
hold on
for i=1:32000
    plot(x(i,1:100),y(i,1:100))
end
hold off


figure
plot(x(:,1),y(:,1))

figure
plot(x(:,50),y(:,50))
figure
plot(x(:,end),y(:,end))


linewidth = 0.5; % desired linewidth 
cmap = jet(19); % colormap, with N colors
x=1:20;
y=1:20;
hold on 
for n = 1:19
    fig=plot(x([n n+1]), y([n n+1]), 'color', cmap(n,:), 'linewidth', linewidth); 
end
colormap(cmap)
colorbar
hold off
