function [] = matcont_plot_label(fig,x,s,dim1,dim2,textFile)
%Label the special points on an eql branch to 1-param bifurcation diagram
% fig: a handle to the figure to plot on
% x, s: output from [x,v,s,h,f]=cont(...)
% dim1: index of the variable for x-axis (ie the bifurcation parameter)
% dim1: index of the variable for y-axis (one of the ODE variables)
% textFile: the text file to save singularity information. 0 for no save
figure(fig);
hold on
if size(s,1)~=2
    % the first and last element of s are start/end point, not interesting
    for i=2:length(s)-1
        type=s(i).label;
        xindex=s(i).index;
        if textFile ~= 0
            fprintf(textFile,['Singularity: type %s at ',repmat('%d, ', 1, size(x,1)),'\n'],type,x(:,xindex));
        end
        fprintf(['Singularity: type %s at ',repmat('%d, ', 1, size(x,1)),'\n'],type,x(:,xindex));
        
        % don't care about neutral saddle (skip fake Hopfs)
        if strcmp(s(i).label,'H ') && strcmp(s(i).msg,'Neutral Saddle Equilibrium')
            continue;
        end
        d=axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3)];
        line(x(dim1,xindex),x(dim2,xindex),'linestyle','none', 'marker', '.','color','r');
        text(x(dim1,xindex)+skew(1),x(dim2,xindex)+skew(2), {s(i).label});
    end
    
end
hold off
end

