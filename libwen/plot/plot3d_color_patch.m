function fig=plot3d_color_patch(x,y,z,v)

hold on
fig=patch(x,y,z,v,'edgecolor','flat','facecolor','none');
view(3);
grid on;
colormap(jet);
colorbar;
hold off

end
