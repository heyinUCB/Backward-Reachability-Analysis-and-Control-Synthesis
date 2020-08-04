% This script is aming at using control Lyapunov function to design
% controller

% discrete time
dt = 0.001;

%% load storage function
load resultT2.mat
DVvalx = jacobian(Vval, x);
DVvalt = jacobian(Vval, t);

Dmax = diag([20 20 50 20]);
% initial condition
x0 = [47; 20; 50; 20];
xeq = [45.0; 0.04924*r2d; 0.0; 0.04924*r2d];
x0 = x0 - xeq;
x0 = inv(Dmax)*x0;

%% Main algorithm
xsim = [];
xsim = [xsim, x0];
uopt = [];
time_line = 0:dt:T;
for i = 1:length(time_line)
    % plug in the measured state and time
    DVx = subs(subs(DVvalx, x, xsim(:,i)),t,time_line(i));
    DVt = subs(subs(DVvalt, x, xsim(:,i)),t,time_line(i));
    
    fval = subs(f,x,xsim(:,i));
    gval = subs(g,x,xsim(:,i));
    a = double(DVx*gval);
    b = double(DVt + DVx*fval);
    % optimal value of u
    if b <= 0 
        uval = 0;
    else
        uval = -b*a'/(a*a');
    end
    uopt = [uopt, uval];
    
    % simulating the system with uval
%     x_next = [xsim(1,i) + dt*xsim(2,i);...
%         
%               xsim(2,i) + dt*(- 10.6560*xsim(1,i)^3 + 11.5309*xsim(1,i)^2*xsim(3,i) + 7.8850*xsim(1,i)*xsim(3,i)^2 + 0.7972*xsim(2,i)^2*xsim(3,i) ...
%               + 0.8408*xsim(2,i)*xsim(3,i)*xsim(4,i) + 21.0492*xsim(3,i)^3 + 0.4204*xsim(3,i)*xsim(4,i)^2 + 66.5225*xsim(1,i) - 24.5110*xsim(3,i)+ ...
%               uval*(-10.0959*xsim(3,i)^2 + 44.2521));...
%               
%               xsim(3,i) + dt*xsim(4,i);...
%               
%               xsim(4,i) + dt*(10.9955*xsim(1,i)^3 - 48.9151*xsim(1,i)^2*xsim(3,i) - 6.4044*xsim(1,i)*xsim(3,i)^2 - 2.3955*xsim(2,i)^2*xsim(3,i) ...
%               - 1.5943*xsim(2,i)*xsim(3,i)*xsim(4,i) - 51.9088*xsim(3,i)^3 - 0.7971*xsim(3,i)*xsim(4,i)^2 - 68.6419*xsim(1,i) + 103.9783*xsim(3,i)+ ...
%               uval*(37.8015*xsim(3,i)^2 - 83.9120))];
    x_next = xsim(:,i)+dt*(fval+gval*uval);
    xsim = [xsim, x_next];
end
%%
xsim = Dmax*xsim;

%%
subplot(2,1,1)
plot(0:dt:T,uopt*100+14.33,'linewidth',2)
xlabel('$t$ (sec)','interpreter','latex')
ylabel('$u \ (\%)$','interpreter','latex')
hold on
grid on

subplot(2,1,2)
plot(0:dt:T,xsim(1,1:end-1)+xeq(1),'linewidth',2)
hold on
plot(0:dt:T,xsim(2,1:end-1)+xeq(2),'linewidth',2)
hold on
plot(0:dt:T,xsim(3,1:end-1)+xeq(3),'linewidth',2)
hold on
plot(0:dt:T,xsim(4,1:end-1)+xeq(4),'linewidth',2)
leg = legend('$U$ (m/s)','$\alpha$ (deg)',...
    '$q$ (deg/s)','$\theta$ (rad/s)');
set(leg,'interpreter','latex')

xlabel('$t$ (sec)','interpreter','latex')
ylabel('states value','interpreter','latex')
grid on
