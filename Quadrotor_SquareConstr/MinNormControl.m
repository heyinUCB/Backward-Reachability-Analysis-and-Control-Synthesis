% This script is aming at using control Lyapunov function to design
% controller

%% System dynamics
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);

x = [x1; x2; x3; x4];
% create time t
t = sdpvar(1,1);

% Define the system
f3 = 0.11707*x2^5 + 0.035908*x2^3*x4^2 - 1.6032*x2^3 - 0.17201*x2*x4^2 + 3.0313*x2;
f4 = 0.24902*x2^5 + 0.13049*x2^3*x4^2 - 5.6188*x2^3 - 0.29147*x2*x4^2 + 23.9892*x2;
f = [x3; x4; f3; f4];
g3 = 0.02905*x2^4 - 0.11289*x2^2 + 0.3955;
g4 = 0.096371*x2^4 - 0.54277*x2^2 + 0.7831;
g = [0; 0; g3; g4];


% define control input
u = sdpvar(1,1);

% discrete time
dt = 0.05;

% time horizon
T = 1;

% initial condition
x0 = [0.15; -0.23; -0.1; 1.06];

%% load storage function
load('Reach15-Aug-2018- 8-07-PM','Vval')
[Vcoeff,R,e] = poly2basis(Vval,Vmonom);

Vdeg = 4;
Vmonom_yalmip = monolist([x; t], Vdeg);
Vval_yalmip = Vcoeff'*Vmonom_yalmip;
DVvalx = jacobian(Vval_yalmip, x);
DVvalt = jacobian(Vval_yalmip, t);

%% Main algorithm
xsim = [];
xsim = [xsim, x0];
uopt = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % plug in the measured state and time
    DVx = replace(replace(DVvalx, x, xsim(:,i)),t,time_line(i));
    DVt = replace(replace(DVvalt, x, xsim(:,i)),t,time_line(i));
    
    options = sdpsettings('verbose',1,'solver','mosek');
    fval = replace(f,x,xsim(:,i));
    gval = replace(g,x,xsim(:,i));
    constr = DVt + DVx*(fval + gval*u) <= 0; 
    obj = u'*u;
    optimize(constr,obj,options)
    % optimal value of u
    uval = value(u);
    uopt = [uopt, uval];
    
    % simulating the system with uval
    x_next = [xsim(1,i) + dt*xsim(3,i);...
              xsim(2,i) + dt*xsim(4,i);...
              xsim(3,i) + dt*(0.11707*xsim(2,i)^5 + 0.035908*xsim(2,i)^3*xsim(4,i)^2 - 1.6032*xsim(2,i)^3 ...
              - 0.17201*xsim(2,i)*xsim(4,i)^2 + 3.0313*xsim(2,i) + uval*(0.02905*xsim(2,i)^4 - 0.11289*xsim(2,i)^2 + 0.3955));...
              xsim(4,i) + dt*(0.24902*xsim(2,i)^5 + 0.13049*xsim(2,i)^3*xsim(4,i)^2 - 5.6188*xsim(2,i)^3 ...
              - 0.29147*xsim(2,i)*xsim(4,i)^2 + 23.9892*xsim(2,i) + uval*(0.096371*xsim(2,i)^4 - 0.54277*xsim(2,i)^2 + 0.7831))];
          
    xsim = [xsim, x_next];
end

%%
figure(1)
plot(0:dt:T,uopt,'linewidth',2)
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u$ (N)','interpreter','latex')

figure(2)
plot(0:dt:T,xsim(1,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(2,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(3,2:end),'linewidth',2)
hold on
plot(0:dt:T,xsim(4,2:end),'linewidth',2)
legend('d (m)','\theta (rad)','v (m/s)','thetadot')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('states value','interpreter','latex')

