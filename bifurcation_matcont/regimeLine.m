function [] = regimeLine(fig,x,s,ind)
%Draws vertical line to separate regimes on 1-parameter bifurcation
%diagrams at a particular bifurcation points
% fig: handle to the bif plot figure
% x, s: output from [x,v,s,h,f]=cont(...)
% ind: the index of the interesting bifurcation point in s

figure(fig);
hold on
xx=x(end,s(ind).index);
line([xx,xx],ylim,'Color',[0.5,0.3,0.5]);
hold off
end

