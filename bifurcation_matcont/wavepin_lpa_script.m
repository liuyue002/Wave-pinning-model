% Bifurcation analysis for wave pinning LPA system
% Using gamma as bifurcation parameter

%% Setup
clear all;clc;

global cds sys
addpath(genpath('matcontfolder'));
addpath(genpath('systemsfolder'));
sys.gui.pausespecial=0;  %Pause at special points 
sys.gui.pausenever=1;    %Pause never 
sys.gui.pauseeachpoint=0; %Pause at each point

% Specify system file, obtain this through the graphical interface
syshandle=@wavepin_wellmix_lpa;

%% Set parameters
SubFunHandles=feval(syshandle);  %Get function handles from system file
RHShandle=SubFunHandles{2};      %Get function handle for ODE

%Set ODE parameter
% gamma is the bifurcation parameter
LL=1^2; w=4.3; k=1*LL; gamma=0.01*LL; n=8; eta=5*LL;

pvec=[w,k,gamma,n,eta]; % Initialize parameter vector
ap=3; % active param
%Set ODE initial condition
uinit=1;
uLinit=1;

RHS_no_param=@(t,x)RHShandle(t,x,w,k,gamma,n,eta); 

%% Find an equilibrium with ODE integration
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);
[tout,xout]=ode45(RHS_no_param,[0,100],[uinit,uLinit],options);
% make sure you have reached an equilibrium. Check by plotting:
%figure;
%plot(tout,xout,'-k','linewidth',2);

%% Continuation on EP branch
% set initial point to be the endpoint of integration (so it's an EP)
xinit=xout(length(xout),:)';

[x0,v0]=init_EP_EP(syshandle, xinit, pvec, ap); %Initialize equilibrium
opt=contset;
opt=contset(opt,'MaxNumPoints',1000); %Set numeber of continuation steps
opt=contset(opt,'MaxStepsize',0.01);  %Set max step size
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues 
opt=contset(opt,'InitStepsize',0.001); %Set Initial stepsize
%opt=contset(opt,'Backward',1);
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); %Equilibrium continuation
% last row of x1 is active param values

%% follow first BP
pointLabel=2; % the label of the branch point in s1
x=x1(1:end-1,s1(pointLabel).index); % u and uL are 1st and 2nd row
pvec(ap)=x1(end,s1(pointLabel).index); % the bif parameter is 3rd row
[x0,v0]=init_BP_EP(syshandle, x, pvec, s1(pointLabel), 0.01);
opt=contset;
opt=contset(opt,'MaxNumPoints',1000); %Set numeber of continuation steps
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'MinStepsize',0.00001);
opt=contset(opt,'Eigenvalues',1); 
opt=contset(opt,'Singularities',1);
%opt=contset(opt,'Multipliers',1); % only for limit cycles
%cds.options.MaxCorrIters=100;
%cds.options.MaxNewtonIters=10;
[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt);

%% follow first BP backwards
% sometimes need to set step size small enough so it follows the right
% branch
pointLabel=2; % the label of the branch point in s1
x=x1(1:end-1,s1(pointLabel).index); % u and uL are 1st and 2nd row
pvec(ap)=x1(end,s1(pointLabel).index); % the bif parameter is 3rd row
[x0,v0]=init_BP_EP(syshandle, x, pvec, s1(pointLabel), 0.01);
opt=contset;
opt=contset(opt,'MaxNumPoints',4000); %Set numeber of continuation steps
opt=contset(opt,'InitStepsize',0.0001);
opt=contset(opt,'MaxStepsize',0.001);
opt=contset(opt,'MinStepsize',0.00001);
opt=contset(opt,'Eigenvalues',1); 
opt=contset(opt,'Backward',1);
opt=contset(opt,'Singularities',1);
%opt=contset(opt,'Multipliers',1); % only for limit cycles
%cds.options.MaxCorrIters=100;
%cds.options.MaxNewtonIters=10;
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt);

%% plotting
titlestr=strcat('wp\_LPA,w=',num2str(w),',k=',num2str(k),',\gamma=',num2str(gamma),',n=',num2str(n),',eta=',num2str(eta));
filename=strrep(titlestr,'\','');
textFile=fopen(strcat(filename,'.txt'),'w');

fig=figure();
matcont_plot_ep(fig,x1,f1,3,2); % 2 is the variable to plot (uL)
matcont_plot_ep(fig,x2,f2,3,2);
matcont_plot_ep(fig,x3,f3,3,2);
xlim([0,3]);
ylim([0,3]);
xlabel('\gamma');
ylabel('u_L');
matcont_plot_label(fig,x1,s1,3,2,textFile);
matcont_plot_label(fig,x2,s2,3,2,textFile);
matcont_plot_label(fig,x3,s3,3,2,textFile);
title(titlestr);
fclose(textFile);

%regimeLineAll(figcopy,x1,s1);
%regimeLineAll(figcopy,x2,s2);
%regimeLineAll(figcopy,x3,s3);

%% Save
saveas(fig,strcat(filename,'.fig'));
%saveas(fig,strcat(filename,'.png'));

