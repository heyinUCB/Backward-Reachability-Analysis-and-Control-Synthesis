%% Reachability Analysis
function [gamma_list, Vval_list, uval] = ...
    CLF(dynamics, T, x, t, rT, V0val, max_iter, uM, um)
warning off

%% load the scaled system
f = dynamics.f;
g = dynamics.g;

%% Define storage function V(t, x) 
Vdeg = 4;
Vmonom = monomials([x; t], 0:Vdeg);
V = polydecvar('V', Vmonom, 'vec');
DVx = jacobian(V, x);
DVt = jacobian(V, t);

%% Define multipliers
sdeg = 4;
% monomials
sxmonom = monomials(x, 0:sdeg);
sxtmonom = monomials([x; t], 0:sdeg);
% multipliers
s1 = polydecvar('s1', sxmonom, 'vec');
s2 = polydecvar('s2', sxtmonom, 'vec');
s3 = polydecvar('s3', sxtmonom, 'vec');
s4 = polydecvar('s4', sxmonom, 'vec');
s5 = polydecvar('s5', sxtmonom, 'vec');
s6 = polydecvar('s6', sxtmonom, 'vec');
s7 = polydecvar('s7', sxtmonom, 'vec');
s8 = polydecvar('s8', sxtmonom, 'vec');
u = polydecvar('u', sxtmonom, 'vec');

%% Define time interval polynomial
w = t*(T-t);

%% Initialization
gamma_list = [];
Vval_list = [];

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
            [s1;s2;s4-0.0001;s6;s8;...
            % V^*(t0,x)<=gamma implies V(t0,x)<= gamma
            -(subs(V,t,0)-gamma_use) + s1*(subs(Vval,t,0)-gamma_use);
            % dV/dt <= 0
            -(DVt + DVx*f + uval*DVx*g) - s2*w + s3val*(V - gamma_use);...
            % V(T,x)<=gamma implies rT<=0
            -s4*rT + subs(V,t,T) - gamma_use;...
            % V(t,x) <= gamma implies u <= uM
            uM - uval + s5val*(V - gamma_use) - s6*w; ...
            % V(t,x) <= gamma implies u >= uM
            uval - um + s7val*(V - gamma_use) - s8*w];
        options = sosoptions;
        options.solver = 'mosek';
        options.solveropts.param.MSK_IPAR_LOG = 0;
        [info,dopt,sossol] = sosopt(SOScons1,[x;t],options);
        if info.feas
            disp('V step is feasible')
            % substitute final decision variable values
            Vval = subs(V,dopt);
            Vval_list = [Vval_list; Vval];
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
    gamma_ub = 1;
    gamma_lb = 0;
    num_exp = 0;
    go = 1;
    while (num_exp <= 12 || go)
        num_exp = num_exp + 1;
        if num_exp >= 20
            break
        end
        gamma_try = (gamma_ub + gamma_lb)/2
        SOScons2 = ...
            [s2;s3;s4-0.0001;s5;s6;s7;s8;...
            % dV/dt <= 0
            -(DVtval + DVxval*f + u*DVxval*g) - s2*w + s3*(Vval - gamma_try); ...
            % V(T,x)<=gamma implies rT<=0
            -s4*rT + subs(Vval,t,T) - gamma_try;...
            % V(t,x) <= gamma implies u <= uM
            uM - u + s5*(Vval - gamma_try) - s6*w; ...
            % V(t,x) <= gamma implies u >= uM
            u - um + s7*(Vval - gamma_try) - s8*w];
        options = sosoptions;
        options.solver = 'mosek';
        options.solveropts.param.MSK_IPAR_LOG = 0;
        [info,dopt,sossol] = sosopt(SOScons2,[x;t],options);
        if info.feas
            disp('g step is feasible')
            gamma_lb = gamma_try;
            go = 0;
            gamma_use = gamma_try;
            % substitute final decision variable values
            uval = subs(u,dopt);
            s3val = subs(s3,dopt);
            s5val = subs(s5,dopt);
            s7val = subs(s7,dopt);
            
        else % if infeasible
            disp('g step is infeasible')
            gamma_ub = gamma_try;
        end
    end
    gamma_list = [gamma_list; gamma_use];
    disp('gamma step ends')
end

