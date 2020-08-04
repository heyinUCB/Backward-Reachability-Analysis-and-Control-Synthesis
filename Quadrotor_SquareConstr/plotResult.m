plotaxes1 = [-1 1 -1.2 1.2];
plotaxes2 = [-1.9 1.9 -1 1];
plotaxes3 = [-pi/10 pi/10 -pi/1.8 pi/1.8];
plotaxes4 = [-2 2 -1 1];

%%
load('V4u2.mat')
gamma_use = gamma_list(end);
V0val = subs(Vval,t,0);

subplot(2,2,1)
x3x4x5x6 = [0;0;0;0];
V012 = subs(V0val,[x3;x4;x5;x6],x3x4x5x6);
[C,h] = pcontour(V012, gamma_use, plotaxes4,'r',[500, 500]);
h.LineColor = mycolor('maroon');
h.LineWidth = 4;
hold on

% VT12 = subs(VTval,[x3;x4;x5;x6],x3x4x5x6);
% pcontour(VT12, gamma_use, plotaxes4,'k',[500, 500]);
% hold on

% rT12 = subs(rt,[x3;x4;x5;x6],x3x4x5x6);
% pcontour(rT12, 0, plotaxes4,'b',[500, 500]);
% hold on

xlabel('$p_x$','interpreter','latex')
ylabel('$p_y$','interpreter','latex')
grid on

%
subplot(2,2,2)
x2x4x5x6 = [0;0;0;0];
V013 = subs(V0val,[x2;x4;x5;x6],x2x4x5x6);
[C,h] = pcontour(V013, gamma_use, plotaxes2,'r',[500, 500]);
h.LineColor = mycolor('maroon');
h.LineWidth = 4;
hold on

% VT13 = subs(VTval,[x2;x4;x5;x6],x2x4x5x6);
% pcontour(VT13, gamma_use, plotaxes2,'k',[500, 500]);
% hold on

% rT13 = subs(rt,[x2;x4;x5;x6],x2x4x5x6);
% pcontour(rT13, 0, plotaxes2,'b',[500, 500]);
% hold on

xlabel('$p_x$','interpreter','latex')
ylabel('$v_x$','interpreter','latex')
grid on

%
subplot(2,2,3)
x1x3x5x6 = [0;0;0;0];
V024 = subs(V0val,[x1;x3;x5;x6],x1x3x5x6);
[C,h] = pcontour(V024, gamma_use, plotaxes1,'r',[500, 500]);
h.LineColor = mycolor('maroon');
h.LineWidth = 4;
hold on

% VT24 = subs(VTval,[x1;x3;x5;x6],x1x3x5x6);
% pcontour(VT24, gamma_use, plotaxes1,'k',[500, 500]);
% hold on

% rT24 = subs(rt,[x1;x3;x5;x6],x1x3x5x6);
% pcontour(rT24, 0, plotaxes1,'b',[500, 500]);
% hold on

xlabel('$p_y$','interpreter','latex')
ylabel('$v_y$','interpreter','latex')
grid on

%
subplot(2,2,4)
x1x2x3x4 = [0;0;0;0];
V056 = subs(V0val,[x1;x2;x3;x4],x1x2x3x4);
[C,h] = pcontour(V056, gamma_use, plotaxes3,'r',[500, 500]);
h.LineColor = mycolor('maroon');
h.LineWidth = 4;
hold on

% VT56 = subs(VTval,[x1;x2;x3;x4],x1x2x3x4);
% pcontour(VT56, gamma_use, plotaxes3,'k',[500, 500]);
% hold on

% rT56 = subs(rt,[x1;x2;x3;x4],x1x2x3x4);
% pcontour(rT56, 0, plotaxes3,'b',[500, 500]);
% hold on

xlabel('$\theta$','interpreter','latex')
ylabel('$\omega$','interpreter','latex')
grid on

%%
load('V2u2.mat')
gamma_use = gamma_list(end);
V0val = subs(Vval,t,0);

subplot(2,2,1)
x3x4x5x6 = [0;0;0;0];
V012 = subs(V0val,[x3;x4;x5;x6],x3x4x5x6);
[C,h] = pcontour(V012, gamma_use, plotaxes4,'r-.',[500, 500]);
h.LineColor = mycolor('orange');
h.LineWidth = 4;
hold on

% VT12 = subs(VTval,[x3;x4;x5;x6],x3x4x5x6);
% pcontour(VT12, gamma_use, plotaxes4,'k',[500, 500]);
% hold on

% rT12 = subs(rt,[x3;x4;x5;x6],x3x4x5x6);
% pcontour(rT12, 0, plotaxes4,'b',[500, 500]);
% hold on

xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
grid on

%
subplot(2,2,2)
x2x4x5x6 = [0;0;0;0];
V013 = subs(V0val,[x2;x4;x5;x6],x2x4x5x6);
[C,h] = pcontour(V013, gamma_use, plotaxes2,'r-.',[500, 500]);
h.LineColor = mycolor('orange');
h.LineWidth = 4;
hold on

% VT13 = subs(VTval,[x2;x4;x5;x6],x2x4x5x6);
% pcontour(VT13, gamma_use, plotaxes2,'k',[500, 500]);
% hold on

% rT13 = subs(rt,[x2;x4;x5;x6],x2x4x5x6);
% pcontour(rT13, 0, plotaxes2,'b',[500, 500]);
% hold on

xlabel('$x_1$','interpreter','latex')
ylabel('$x_3$','interpreter','latex')
grid on

%
subplot(2,2,3)
x1x3x5x6 = [0;0;0;0];
V024 = subs(V0val,[x1;x3;x5;x6],x1x3x5x6);
[C,h] = pcontour(V024, gamma_use, plotaxes1,'r-.',[500, 500]);
h.LineColor = mycolor('orange');
h.LineWidth = 4;
hold on

% VT24 = subs(VTval,[x1;x3;x5;x6],x1x3x5x6);
% pcontour(VT24, gamma_use, plotaxes1,'k',[500, 500]);
% hold on

% rT24 = subs(rt,[x1;x3;x5;x6],x1x3x5x6);
% pcontour(rT24, 0, plotaxes1,'b',[500, 500]);
% hold on

xlabel('$x_2$','interpreter','latex')
ylabel('$x_4$','interpreter','latex')
grid on

%
subplot(2,2,4)
x1x2x3x4 = [0;0;0;0];
V056 = subs(V0val,[x1;x2;x3;x4],x1x2x3x4);
[C,h] = pcontour(V056, gamma_use, plotaxes3,'r-.',[500, 500]);
h.LineColor = mycolor('orange');
h.LineWidth = 4;
hold on

% VT56 = subs(VTval,[x1;x2;x3;x4],x1x2x3x4);
% pcontour(VT56, gamma_use, plotaxes3,'k',[500, 500]);
% hold on

% rT56 = subs(rt,[x1;x2;x3;x4],x1x2x3x4);
% pcontour(rT56, 0, plotaxes3,'b',[500, 500]);
% hold on

xlabel('$x_5$','interpreter','latex')
ylabel('$x_6$','interpreter','latex')
grid on

%%
subplot(2,2,1)
box1 = Polyhedron('ub',[1.7; 0.85], 'lb',[-1.7; -0.85]);
box1.plot('alpha',0.1,'color','b','linewith',4)
legend('deg 4', 'deg 2', '$\Omega_{t,0}^r$','interpreter','latex')

subplot(2,2,2)
box2 = Polyhedron('ub',[1.7; 0.8], 'lb',[-1.7; -0.8]);
box2.plot('alpha',0.1,'color','b','linewith',4)

subplot(2,2,3)
box3 = Polyhedron('ub',[0.85; 1], 'lb',[-0.85; -1]);
box3.plot('alpha',0.1,'color','b','linewith',4)

subplot(2,2,4)
box4 = Polyhedron('ub',[pi/12; pi/2], 'lb',[-pi/12; -pi/2]);
box4.plot('alpha',0.1,'color','b','linewith',4)

garyfyFigure
