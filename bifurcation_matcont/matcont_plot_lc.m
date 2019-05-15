function [] = matcont_plot_lc(fig,x,f,dim1,dim2,num_dim)
%Plot an limit cycle branch from matcont on 1-param bifurcation diagram
% fig: a handle to the figure to plot on
% x, f: output from [x,v,s,h,f]=cont(...)
% Make sure set 'Multipliers' to 1, so f will contain Floquet multipliers
% for determining stability of limit cycles
% dim1: index of the variable for x-axis (ie the bifurcation parameter)
% dim2: index of the variable for y-axis (one of the ODE variables)
% num_dim: number of variables in the original ODE (do not count bif param)
figure(fig);
hold on
% the bif param is always the last row of x
% the 2nd last row are periods
lcmax=max(x(dim2:num_dim:end-2,:));
lcmin=min(x(dim2:num_dim:end-2,:));
%plot(x2(end,:),lcmax);
%plot(x2(end,:),lcmin);
for i=1:size(x,2)
    if abs(f(end-1,i)*f(end,i))<1 % the Lyapunov exponent
        %limit cycle is stable
        plot(x(end,i),lcmax(i),'og');
        plot(x(end,i),lcmin(i),'og');
    else
        %limit cycle is unstable
        plot(x(end,i),lcmax(i),'ob');
        plot(x(end,i),lcmin(i),'ob');
    end
end

hold off
end

