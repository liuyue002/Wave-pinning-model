% Bifurcation analysis for non-conservative LPA system
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
syshandle=@champ_lpa;

%% Set parameters
SubFunHandles=feval(syshandle);  %Get function handles from system file
RHShandle=SubFunHandles{2};      %Get function handle for ODE

%Set ODE parameter
% gamma is the bifurcation parameter
LL=1^2; k=1*LL; gamma=0.01*LL; n=2; eta=5*LL; epc=1; alpha=1.5*LL; theta=5.5*LL;

pvec=[k,gamma,n,eta,epc,alpha,theta]; % Initialize parameter vector
ap=2; %active param
%Set ODE initial condition
uinit=1;
vinit=1;
uLinit=1;

RHS_no_param=@(t,x)RHShandle(t,x,k,gamma,n,eta,epc,alpha,theta); 

%% Find an equilibrium with ODE integration
options=odeset;
options=odeset(options,'RelTol',1e-5);
options=odeset(options,'maxstep',1e-1);
[tout,xout]=ode45(RHS_no_param,[0,1000],[uinit,vinit,uLinit],options);
% make sure you have reached an equilibrium. Check by plotting:
%figure;
%plot(tout,xout,'-k','linewidth',2);

%% Continuation on EP branch
% set initial point to be the endpoint of integration (so it's an EP)
xinit=xout(length(xout),:)';

[x0,v0]=init_EP_EP(syshandle, xinit, pvec, ap); %Initialize equilibrium
opt=contset;
opt=contset(opt,'MaxNumPoints',5000); %Set numeber of continuation steps
opt=contset(opt,'MaxStepsize',0.02);  %Set max step size
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues 
opt=contset(opt,'InitStepsize',0.01); %Set Initial stepsize
opt=contset(opt,'Backward',1); % for this one backward is forward
[x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); %Equilibrium continuation
% last row of x1 is active param values

%% Follow branch point forward (get low branch)
% It seems that always should do eqlib branches before limit cycles,
% otherwise can get error
pointLabel=2; % the label of the branch point in s1
x=x1(1:3,s1(pointLabel).index); % u and uL are 1st and 2nd row
pvec(ap)=x1(4,s1(pointLabel).index); % the bif parameter is 4rd row
[x0,v0]=init_BP_EP(syshandle, x, pvec, s1(pointLabel), 0.01);
opt=contset;
opt=contset(opt,'MaxNumPoints',5000); %Set numeber of continuation steps
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MaxStepsize',0.05);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'Eigenvalues',1); 
[x2,v2,s2,h2,f2]=cont(@equilibrium,x0,v0,opt);

%% Follow branch point backward (get high branch)
pointLabel=2; % the label of the branch point in s1
x=x1(1:3,s1(pointLabel).index); % u and uL are 1st and 2nd row
pvec(ap)=x1(4,s1(pointLabel).index); % the bif parameter is 4rd row
[x0,v0]=init_BP_EP(syshandle, x, pvec, s1(pointLabel), 0.001);
opt=contset;
opt=contset(opt,'MaxNumPoints',10000); %Set numeber of continuation steps
opt=contset(opt,'Singularities',1);  %Monitor singularities
opt=contset(opt,'InitStepsize',0.001);
opt=contset(opt,'MaxStepsize',0.05);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'Eigenvalues',1); 
opt=contset(opt,'Backward',1);
[x3,v3,s3,h3,f3]=cont(@equilibrium,x0,v0,opt);

% the first point might be out of place for some reason, can delete it
%x3(:,1)=x3(:,2);

%% Follow main branch Hopf
pointLabel=3; % the label of the Hopf point in s1
xHopf=x1(1:3,s1(pointLabel).index); % u and L are 1st and 2nd row
pvec(ap)=x1(4,s1(pointLabel).index); % the bif parameter is 4rd row

%last 3 params are: h, ntst, ncol (numerics settings)
[x0,v0]=init_H_LC(syshandle, xHopf, pvec, ap, 1e-4, 15, 4);
opt=contset;
opt=contset(opt,'MaxNumPoints',5000);
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'Singularities',1);
[x4,v4,s4,h4,f4]=cont(@limitcycle,x0,v0,opt);

%% Follow low branch Hopf
pointLabel=2; % the label of the Hopf point in s1
xHopf=x2(1:3,s2(pointLabel).index); % u and L are 1st and 2nd row
pvec(ap)=x2(4,s2(pointLabel).index); % the bif parameter is 4rd row

[x0,v0]=init_H_LC(syshandle, xHopf, pvec, ap, 1e-4, 15, 4);
opt=contset;
opt=contset(opt,'MaxNumPoints',5000);
opt=contset(opt,'InitStepsize',0.1);
opt=contset(opt,'MaxStepsize',0.5);
opt=contset(opt,'MinStepsize',0.000001);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'Singularities',1);
[x5,v5,s5,h5,f5]=cont(@limitcycle,x0,v0,opt);

%% Follow high branch Hopf
pointLabel=3; % the label of the Hopf point in s1
xHopf=x3(1:3,s3(pointLabel).index); % u and L are 1st and 2nd row
pvec(ap)=x3(4,s3(pointLabel).index); % the bif parameter is 4rd row

[x0,v0]=init_H_LC(syshandle, xHopf, pvec, ap, 1e-2, 15, 4);
opt=contset;
opt=contset(opt,'MaxNumPoints',5000);
opt=contset(opt,'InitStepsize',0.01);
opt=contset(opt,'MaxStepsize',0.1);
opt=contset(opt,'MinStepsize',0.001);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'Singularities',1);
[x6,v6,s6,h6,f6]=cont(@limitcycle,x0,v0,opt);

%% plotting
titlestr=strcat('champ\_LPA',',k=',num2str(k),',\gamma=',num2str(gamma),',n=',num2str(n),',\eta=',num2str(eta),',\epsilon_c=',num2str(epc),',\alpha=',num2str(alpha),',\theta=',num2str(theta));
filename=strrep(titlestr,'\','');
textFile=fopen(strcat(filename,'.txt'),'w');

fig=figure();
xlim([0,100]); ylim([0,3.5]);
matcont_plot_ep(fig,x1,f1,4,3);
matcont_plot_ep(fig,x2,f2,4,3);
matcont_plot_ep(fig,x3,f3,4,3);
matcont_plot_lc(fig,x4,f4,4,3,3);
matcont_plot_lc(fig,x5,f5,4,3,3);
matcont_plot_lc(fig,x6,f6,4,3,3);
xlabel('\gamma');
ylabel('u_L');
matcont_plot_label(fig,x1,s1,4,3,textFile);
matcont_plot_label(fig,x2,s2,4,3,textFile);
matcont_plot_label(fig,x3,s3,4,3,textFile);
title(titlestr);
fclose(textFile);

%regimeLineAll(figcopy,x1,s1);
%regimeLineAll(figcopy,x2,s2);
%regimeLineAll(figcopy,x3,s3);

%% Save
saveas(fig,strcat(filename,'.fig'));

