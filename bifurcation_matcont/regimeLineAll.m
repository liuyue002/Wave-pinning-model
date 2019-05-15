function [] = regimeLineAll(fig,x,s)
%Draws vertical line to separate regimes on 1-parameter bifurcation
%diagrams at ALL bifurcation points
% fig: handle to the bif plot figure
% x, s: output from [x,v,s,h,f]=cont(...)

if length(s) <= 2
    return;
end

figure(fig);
hold on
for i=2:length(s)-1
    % don't care about neutral saddle
    if strcmp(s(i).label,'H ') && strcmp(s(i).msg,'Neutral Saddle Equilibrium')
        continue;
    end
    xx=x(end,s(i).index);
    line([xx,xx],ylim,'Color',[0.5,0.3,0.5]);
end
hold off
end

