function [] = matcont_plot_ep(fig,x,f,dim1,dim2)
%Plot an equilibrium branch from matcont on 1-param bifurcation diagram
% fig: a handle to the figure to plot on
% x, f: output from [x,v,s,h,f]=cont(...)
% Make sure set 'Eigenvalues' to 1, so f will contain eigenvalues for
% determining stability of equilibriums
% dim1: index of the variable for x-axis (ie the bifurcation parameter)
% dim2: index of the variable for y-axis (one of the ODE variables)
figure(fig);
hold on
%plot(x(dim1,:),x(dim2,:),'-r');
for i=2:size(x,2)
    if f(:,i)<0 %the eigenvalue
        %ep is stable
        plot(x(dim1,i-1:i),x(dim2,i-1:i),'-r');
    else
        %ep is unstable
        plot(x(dim1,i-1:i),x(dim2,i-1:i),'-k');
    end
end

hold off
end

