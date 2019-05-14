function [fig] = plot_kymograph(uu, fig_pos,T)
%Plot kymograph for PDE simulation results
% uu: matrix where uu(i,j) contains the value at t_i, x_j
% fig_pos: a vector of 4 numbers, specify location and size of figure
% T: total time of simulation
nFrame=size(uu,1);
nx=size(uu,2);
tTick=(0:0.2:1)*nFrame;
tTickLabel=num2cell((0:0.2:1)*T);
tTickLabel=cellfun(@num2str,tTickLabel,'un',0);
xTick=(0:0.2:1)*nx;
xTickLabel=num2cell((0:0.2:1));
xTickLabel=cellfun(@num2str,xTickLabel,'un',0);

fig=figure('Position',fig_pos);
axis([0 nFrame 0 nx]);
umax=max(max(uu(10:end,:)));
umax=ceil(umax*10)/10;
umin=min(min(uu(10:end,:)));
umin=floor(umin*10)/10;
ucolortick=[umin,umax];
imagesc(uu',ucolortick);
set(gca,'YDir','normal');
colorbar('FontSize',40,'TickLabels',ucolortick,'Ticks',ucolortick);
set(gca,'XTick',tTick);
set(gca,'XTickLabel',tTickLabel);
set(gca,'YTick',xTick);
set(gca,'YTickLabel',xTickLabel);
xlabel('t');
ylabel('x');
%biggerFont(gca);
%tightfig(fig);

end

