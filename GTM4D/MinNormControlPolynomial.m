% This script is aming at using control Lyapunov function to design
% controller

% discrete time
dt = 0.002;
r2d = 180/pi;

%% load storage function
load step4.mat;
Vval = Vval_list(end);
DVvalx = jacobian(Vval, x);
DVvalt = jacobian(Vval, t);

Dmax = diag([20 20 50 20]);
% initial condition
x0 = [47; 20; 100; 20];
xeq = [45.0; 0.04924*r2d; 0.0; 0.04924*r2d];
x0 = x0 - xeq;
x0 = inv(Dmax)*x0;

%% Main algorithm
xsim = [];
xsim = [xsim, x0];
uopt = [];
Vtx = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % plug in the measured state and time    
    fval = subs(f,x,xsim(:,i));
    gval = subs(g,x,xsim(:,i));
    Vtx = [Vtx; subs(subs(Vval,x,xsim(:,i)),t,time_line(i))];

    ucmd = subs(subs(uval, x, xsim(:,i)),t,time_line(i));
    uopt = [uopt, ucmd];
    
    % simulating the system with uval
    x_next = xsim(:,i)+dt*(fval+gval*ucmd);     
    xsim = [xsim, x_next];
end

%%
xsim = Dmax*xsim;

%%
subplot(2,1,1)
plot(0:dt:T,(uopt+0.0489)*180/pi,'linewidth',2)
hold on
plot(0:dt:T, 0.0489*180/pi*ones(1,size(0:dt:T,2)),'-.','linewidth',2)
hold on
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u$ (deg)','interpreter','latex')
grid on
leg = legend('$u_{elev}$ (deg)','$u_{elev,t}$ (deg)');
set(leg,'interpreter','latex')

subplot(2,1,2)
plot(0:dt:T,xsim(1,1:end-1)+xeq(1),'linewidth',2)
hold on
plot(0:dt:T,xsim(2,1:end-1)+xeq(2),'linewidth',2)
hold on
plot(0:dt:T,xsim(3,1:end-1)+xeq(3),'linewidth',2)
hold on
plot(0:dt:T,xsim(4,1:end-1)+xeq(4),'linewidth',2)
leg = legend('$x_1$ (m/s)','$x_2$ (deg)',...
    '$x_3$ (deg/s)','$x_4$ (deg)');
set(leg,'interpreter','latex')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('states value','interpreter','latex')
grid on
