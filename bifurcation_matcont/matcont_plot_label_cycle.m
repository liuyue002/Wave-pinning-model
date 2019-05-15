function [] = matcont_plot_label_cycle(fig,x,s,dim1,dim2,num_dim,textFile)
%Label the special points on a cycle to 1-param bifurcation diagram
% fig: a handle to the figure to plot on
% x, s: output from [x,v,s,h,f]=cont(...)
% dim1: index of the variable for x-axis (ie the bifurcation parameter) 
%(ignored, it's always the last row in x)
% dim2: index of the variable for y-axis (one of the ODE variables)
% num_dim: number of variables in the original ODE (do not count bif param)
% textFile: the text file to save singularity information. 0 for no save
figure(fig);
hold on
if size(s,1)~=2
    % the first and last element of s are start/end point, not interesting
    for i=2:length(s)-1
        type=s(i).label;
        xindex=s(i).index;
        lcmax=max(x(dim2:num_dim:end-2,xindex));
        %lcmin=min(x(dim2:num_dim:end-2,xindex));
        if textFile ~= 0
            fprintf(textFile,'Singularity: type %s at var(cyclemax) = %d, period = %d, par=%d\n',type,lcmax,x(end-1:end,xindex));
        end
        fprintf('Singularity: type %s at var(cyclemax) = %d, period = %d, par=%d\n',type,lcmax,x(end-1:end,xindex));
        
        d=axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3)];
        line(x(end,xindex),lcmax,'linestyle','none', 'marker', '.','color','r');
        text(x(end,xindex)+skew(1),lcmax+skew(2), {s(i).label});
    end
    
end
hold off
end

