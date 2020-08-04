function pcontourslice(p,v,domain,npts,var,xslice,yslice,zslice)

%-------------------------------------------------------------------
% Set defaults 
%-------------------------------------------------------------------

% Default contour
if isempty(v)
    v=1;
end
lv = length(v);

% Default npts
if isempty(npts)    
    Nx = 100;
    Ny = 100;
    Nz = 5;
else
    Nx = npts(1);
    Ny = npts(2);
    Nz = npts(3);
end

% Define variables as chars
if isempty(var)
    x = p.varname{1};
    y = p.varname{2};
    z = p.varname{3};
elseif ispvar(var) || iscellstr(var)
    if ispvar(var)
        var = char(var);
    end
    x = var{1};
    y = var{2};
    z = var{3};
else
    error('var must be a vector of pvars');
end

if isempty(domain)
    domain = [-1 1 -1 1 -1 1];
end    
    
% Plot contour
xg = linspace(domain(1),domain(2),Nx);
yg = linspace(domain(3),domain(4),Ny);
zg = linspace(domain(5),domain(6),Nz);
[xg,yg,zg] = meshgrid(xg,yg,zg);
%pgrid = double(subs(p,{x,y},{xg,yg}));
pgrid = double(subs(p,{x; y; z},[xg(:)'; yg(:)'; zg(:)']));
pgrid = reshape(pgrid,size(xg));
if lv==1
    % Single contour syntax for contour function
    v = [v v];
end

contourslice(xg,yg,zg,pgrid,xslice,yslice,zslice,v);    
xlabel(x)
ylabel(y)
zlabel(z)
xlim([domain(1),domain(2)])
ylim([domain(3),domain(4)])
zlim([domain(5),domain(6)])
