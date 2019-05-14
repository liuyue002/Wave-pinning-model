% Finite difference simulation for the extended 
% wave pinning model in one spatial dimension

%% Numerics parameters
drawperframe=50; % update the animation once per 50 iterations
nx=400; % number of spatial discretization points
dx=1/nx;
T=1; % total time for the simulation
dt=0.0002; % time discretization
nt=T/dt+1;
nFrame=ceil((T/dt)/drawperframe);

%% Model parameters
LL=20^2; delta=0.01; k0=1.0*LL; gamma=18.5*LL; n=2; eta=5.2*LL;
c=1; theta=5.5*LL; alpha=1.5*LL;
s=40*LL; epsilon=0.1; kn=24*LL; ks=7.5*LL;

%% Set up

f = @(u,v,F) (k0+gamma.*u.^n./(1 + u.^n)).*v - (eta+s*F./(1 + F)).*u;

% the unique equilibrium
u0=alpha/theta;
f0=u0*kn/ks;
v0=(c*alpha+(eta+s*f0/(1+f0))*u0)/(k0+gamma*u0^n/(1+u0^n));

x=linspace(0,1,nx)';
u=zeros(nx,1); v=zeros(nx,1); F=zeros(nx,1);

o=ones(nx,1);
A=spdiags([o -2*o o],[-1 0 1],nx,nx);
A(1,1)=-1; % for no-flux BC
A(nx,nx)=-1;
A=A/(dx^2); % A is the discrete Laplacian

%% initial conditions
u=u0*((rand(size(u)))*0.2+0.9);
v(:)=v0;
F(:)=f0;

%% Set up figure
fig=figure();
xlabel('x'); ylabel('u,v,F');
axis([0 1 0 6]);

uu=zeros(nFrame,nx); %for plotting kymographs
vv=zeros(nFrame,nx);
ff=zeros(nFrame,nx);

%% Main solver
th=0.5; % 0: forward Euler, 0.5: Crank-Nicolson, 1: backward Euler
Tu=speye(nx)-th*dt*delta*A;
Tv=speye(nx)-th*dt*A;
for ti=1:1:nt
    if (mod(ti, drawperframe) == 1)
        cla;
        hold on;
        plot(x,u); plot(x,v); plot(x,F);
        legend('u','v','F');
        hold off;
        drawnow
        uu(iFrame,:)=u;
        vv(iFrame,:)=v;
        ff(iFrame,:)=F;
    end
    
    urhs = u + dt*(f(u,v,F) + (1-th)*delta*A*u - c*theta*u);
    unew = Tu\urhs;
    
    vrhs = v + dt*(-f(u,v,F) + (1-th)*A*v + c*alpha);
    vnew = Tv\vrhs;
    
    Fnew = F +(dt*epsilon)*(kn*u - ks*F);
    F=Fnew; u=unew; v=vnew;
end

figu=plot_kymograph(uu, fig_pos, T);
figv=plot_kymograph(vv, fig_pos, T);
figF=plot_kymograph(ff, fig_pos, T);
