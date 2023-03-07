%% Reachability Analysis
function [gamma_list, Vval, u1val, u2val] = ...
    CLF(dynamics, T, x, t, rt, V0val, max_iter, uM, um, epsilon)
warning off

%% load the bound of input
uM1 = uM(1);
um1 = um(1);
uM2 = uM(2);
um2 = um(2);

%% load the system
f = dynamics.f;
g1 = dynamics.g1;
g2 = dynamics.g2;

%% load the state constraints
rt1 = rt.rt1;
rt2 = rt.rt2;
rt3 = rt.rt3;
rt4 = rt.rt4;
rt5 = rt.rt5;
rt6 = rt.rt6;

%% Define storage function V(t,x)
Vdeg = 2;
Vmonom = monomials([x; t], 0:Vdeg);
V = polydecvar('V', Vmonom, 'vec');
DVx = jacobian(V, x);
DVt = jacobian(V, t);

%% Define multipliers
sdeg = 2;
% monomials
sxmonom = monomials(x, 0:sdeg);
sxtmonom = monomials([x; t], 0:sdeg);
% multipliers
s1 = polydecvar('s1', sxmonom, 'vec');
s2 = polydecvar('s2', sxtmonom, 'vec');
s3 = polydecvar('s3', sxtmonom, 'vec');
s4a = polydecvar('s4a', sxtmonom, 'vec');
s4b = polydecvar('s4b', sxtmonom, 'vec');
s4c = polydecvar('s4c', sxtmonom, 'vec');
s4d = polydecvar('s4d', sxtmonom, 'vec');
s4e = polydecvar('s4e', sxtmonom, 'vec');
s4f = polydecvar('s4f', sxtmonom, 'vec');
s5a = polydecvar('s5a', sxtmonom, 'vec');
s5b = polydecvar('s5b', sxtmonom, 'vec');
s6a = polydecvar('s6a', sxtmonom, 'vec');
s6b = polydecvar('s6b', sxtmonom, 'vec');
s7a = polydecvar('s7a', sxtmonom, 'vec');
s7b = polydecvar('s7b', sxtmonom, 'vec');
s8a = polydecvar('s8a', sxtmonom, 'vec');
s8b = polydecvar('s8b', sxtmonom, 'vec');
s9a = polydecvar('s9a', sxtmonom, 'vec');
s9b = polydecvar('s9b', sxtmonom, 'vec');
s9c = polydecvar('s9c', sxtmonom, 'vec');
s9d = polydecvar('s9d', sxtmonom, 'vec');
s9e = polydecvar('s9e', sxtmonom, 'vec');
s9f = polydecvar('s9f', sxtmonom, 'vec');
u1 = polydecvar('u1', sxtmonom, 'vec');
u2 = polydecvar('u2', sxtmonom, 'vec');

%% Define time interval polynomial
w = t*(T-t);

%% Initialization
gamma_list = [];

for i = 1:max_iter
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    fprintf(['This is the ',num2str(i),'th iteration\n'])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%% V step %%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('V step starts')
    
    if i == 1
        Vval = V0val;
    else
        % gamma is fixed
        SOScons1 = ...
            [s1;s2;...
            s4a-epsilon;...
            s4b-epsilon;...
            s4c-epsilon;...
            s4e-epsilon;...
            s4d-epsilon;...
            s4e-epsilon;...
            s6a;s6b;s8a;s8b;s9a;s9b;s9c;s9d;s9e;s9f;...
            % V^*(t0,x)<=gamma implies V(t0,x)<= gamma
            -(subs(V,t,0)-gamma_use) + s1*(subs(Vval,t,0) - gamma_use);
            % dV/dt <= 0
            -(DVt + DVx*f + u1val*DVx*g1 + u2val*DVx*g2) - s2*w ...
                                        + s3val*(V - gamma_use);...
            % V(t,x)<=gamma implies rt<=0
            s4a*rt1 + (V - gamma_use) - s9a*w;...
            s4b*rt2 + (V - gamma_use) - s9b*w;...
            s4c*rt3 + (V - gamma_use) - s9c*w;...
            s4d*rt4 + (V - gamma_use) - s9d*w;...
            s4e*rt5 + (V - gamma_use) - s9e*w;...
            s4f*rt6 + (V - gamma_use) - s9f*w;...
            % V(t,x) <= gamma implies u <= uM
            uM1 - u1val + s5aval*(V - gamma_use) - s6a*w;...
            uM2 - u2val + s5bval*(V - gamma_use) - s6b*w;...
            % V(t,x) <= gamma implies u >= uM
            u1val - um1 + s7aval*(V - gamma_use) - s8a*w;...
            u2val - um2 + s7bval*(V - gamma_use) - s8b*w];
        options = sosoptions;
        options.solver = 'mosek';
        options.solveropts.param.MSK_IPAR_LOG = 0;
        [info,dopt,sossol] = sosopt(SOScons1,[x;t],options);
        if info.feas
            disp('V step is feasible')
            % substitute final decision variable values
            Vval = subs(V,dopt);
        else
            disp('V step is infeasible')
            break
        end
    end
    disp('V step ends')
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%% gamma step %%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('gamma step starts')
    % V is given, maximize gamma 
    % Solve the SOS feasibility problem
    DVxval = jacobian(Vval, x);
    DVtval = jacobian(Vval, t);
    
    % Initialize the bisection
    gamma_ub = 4;
    gamma_lb = 0;
    num_exp = 0;
    go = 1;
    while (num_exp <= 12 || go)
        num_exp = num_exp + 1;
        gamma_try = (gamma_ub + gamma_lb)/2
        SOScons2 = ...
            [s2;s3;...
            s4a-epsilon;...
            s4b-epsilon;...
            s4c-epsilon;...
            s4e-epsilon;...
            s4d-epsilon;...
            s4e-epsilon;...
            s5a;s5b;s6a;s6b;s7a;s7b;s8a;s8b;...
            s9a;s9b;s9c;s9d;s9e;s9f;...
            % dV/dt <= 0
            -(DVtval + DVxval*f + u1*DVxval*g1 + u2*DVxval*g2) - s2*w ...
                                              + s3*(Vval - gamma_try); ...
            % V(t,x)<=gamma implies rt<=0
            s4a*rt1 + (Vval - gamma_try) - s9a*w;...
            s4b*rt2 + (Vval - gamma_try) - s9b*w;...
            s4c*rt3 + (Vval - gamma_try) - s9c*w;...
            s4d*rt4 + (Vval - gamma_try) - s9d*w;...
            s4e*rt5 + (Vval - gamma_try) - s9e*w;...
            s4f*rt6 + (Vval - gamma_try) - s9f*w;...
            % V(t,x) <= gamma implies u <= uM
            uM1 - u1 + s5a*(Vval - gamma_try) - s6a*w; ...
            uM2 - u2 + s5b*(Vval - gamma_try) - s6b*w; ...
            % V(t,x) <= gamma implies u >= um
            u1 - um1 + s7a*(Vval - gamma_try) - s8a*w;...
            u2 - um2 + s7b*(Vval - gamma_try) - s8b*w];
        options = sosoptions;
        options.solver = 'mosek';
        options.solveropts.param.MSK_IPAR_LOG = 0;
        [info,dopt,sossol] = sosopt(SOScons2,[x;t],options);
        if info.feas
            disp('gamma step is feasible')
            gamma_lb = gamma_try;
            go = 0;
            gamma_use = gamma_try;
            % substitute final decision variable values
            u1val = subs(u1,dopt);
            u2val = subs(u2,dopt);
            s3val = subs(s3,dopt);
            s5aval = subs(s5a,dopt);
            s5bval = subs(s5b,dopt);
            s7aval = subs(s7a,dopt);
            s7bval = subs(s7b,dopt);
            
        else % if infeasible
            disp('gamma step is infeasible')
            gamma_ub = gamma_try;
        end
    end
    gamma_list = [gamma_list; gamma_use];
    disp('gamma step ends')
end

