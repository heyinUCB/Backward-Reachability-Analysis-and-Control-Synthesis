clear
format long
 
%%
x = sdpvar(6,1);
x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4); 
x5 = x(5); 
x6 = x(6); 
g = 9.81;
L = 0.25;
I = 0.00383;
m = 0.486;
f = [x2; 0; x4; -g; x6; 0];
g1 = [0; (0.16*x5^3-x5)/m; 0; (-0.47*x5^2+1)/m; 0; -L/I];
g2 = [0; (0.16*x5^3-x5)/m; 0; (-0.47*x5^2+1)/m; 0; L/I];

%% Define storage function V(t, x) = V0(x) + tV1(x)
Vdeg = 2;
% Vmonom = monolist(x, Vdeg);
% Vcoeff = sdpvar(size(Vmonom,1),1);
% V = Vcoeff'*Vmonom;
[V, Vcoeff] = polynomial(x, Vdeg);
DVx = jacobian(V, x);

%% Define s(x) and s'(x)
sdeg = 2;
[l1,cl1] = polynomial(x, sdeg);
[l2,cl2] = polynomial(x, sdeg);
[s1, c1] = polynomial(x, sdeg);

%% Main algorithm
SOScons = [sos(-DVx*f + l1*DVx*g1 + l2*DVx*g2);...
    sos(V - 0.0001)...
    ];

obj = [];
options = sdpsettings('verbose',1,'solver','penbmi');
options.penbmi.PBM_MAX_ITER = 50;
solvesos(SOScons,obj,options,[Vcoeff;cl1;cl2])

%%
V0coeffval = value(Vcoeff);