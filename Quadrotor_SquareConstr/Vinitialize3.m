format long
 
%%
x = sdpvar(6,1);
x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4); 
x5 = x(5); 
x6 = x(6); 
t = sdpvar(1,1);
alpha_dec = sdpvar(1,1);

% parameters
g = 9.81;
L = 0.25;
I = 0.00383;
m = 0.486;
f = [x2; 0; x4; -g; x6; 0];
g1 = [0; (0.16*x5^3-x5)/m; 0; (-0.47*x5^2+1)/m; 0; -L/I];
g2 = [0; (0.16*x5^3-x5)/m; 0; (-0.47*x5^2+1)/m; 0; L/I];

%% Define storage function
Vdeg = 2;
% Vmonom = monolist(x, Vdeg);
% Vcoeff = sdpvar(size(Vmonom,1),1);
% V = Vcoeff'*Vmonom;
xbar = [x1^2;x1*x2;x2^2;x1*x3;x2*x3;x3^2;x1*x4;x2*x4;x3*x4;x4^2;x1*x5;...
    x2*x5;x3*x5;x4*x5;x5^2;x1*x6;x2*x6;x3*x6;x4*x6;x5*x6;x6^2];
M = sdpvar(1,size(xbar,1));
V = x'*P*x + xbar'*diag(M)*xbar;

DVx = jacobian(V, x);

%% Define s(x) and s'(x)
sdeg = 2;
[l1,cl1] = polynomial(x, sdeg);
[l2,cl2] = polynomial(x, sdeg);
[s1,c1] = polynomial(x, sdeg);
[s2,c2] = polynomial(x, sdeg);
[s3,c3] = polynomial(x, sdeg);
[s4,c4] = polynomial(x, sdeg);

%% local polynomials
T = 1;
w = t*(T-t);

%% Main algorithm
% SOScons = [sos(-(s1*DVx*f + 0.0001 + l*DVx*g))...
%     sos(V - 0.0001);...
%     sos(s1)];
% 
% obj = M*M';
% options = sdpsettings('verbose',1,'solver','penbmi');
% solvesos(SOScons,obj,options,[M';cl;c1])

SOScons = [sos(-DVx*f + l1*DVx*g1 + l2*DVx*g2);...
    sos(V - 0.0001)];

obj = M*M';
options = sdpsettings('verbose',1,'solver','penbmi');
options.penbmi.PBM_MAX_ITER = 50;
solvesos(SOScons,obj,options,[M';cl1;cl2;c1])

%% get the optimal value of decision vars
Mval = value(M);
save Mval.mat Mval