function moment_vec = MomentOfMeasure(nX, qdeg, ubvec, lbvec)
% nX = 6;
% qdeg = 6;
% ubvec = [1.5;1.5;1.5;1.5;1.5;1.5];
% lbvec = [-1.5;-1.5;-1.5;-1.5;-1.5;-1.5];
pvar x1 x2 x3 x4 x5 x6
statevec = [x1; x2; x3; x4; x5; x6];
% choose the states according to the number of states
x = statevec(1:nX);
full_index_list = 1:nX;

% q(x) is the polynomial used to approximate the indicator function
qmonom = monomials(x, 0:qdeg);

% number of monomials of q(x)
num_monom = size(qmonom,1);

% initialize the momoent of measure vector
moment_vec = [];

%% loop through all the monomials
for k = 1:num_monom
    a = qmonom(k,1);
    
    % number of polynomial variables in each monomials
    num_variable = length(a.varname);
    
    % initialize an index list
    index_list = [];
    % initialize one moment
    moment = 1;
    for j = 1:num_variable
        for i = 1:nX
            if contains(a.varname(j),num2str(i))
                % find the variable's index in the state vector
                index = i;
                index_list = [index_list; index];
                % find its degree
                deg = full(a.degmat(j));
                moment = moment*...
                    1/(deg+1)*(ubvec(index)^(deg+1) - lbvec(index)^(deg+1));
                continue
            end
        end
    end
    complement_index_list = setdiff(full_index_list, index_list);
    for m = 1:length(complement_index_list)
        complement_index = complement_index_list(m);
        moment = moment*(ubvec(complement_index)-lbvec(complement_index));
    end
    moment_vec = [moment_vec; moment];
end
