function Vcoeff = loadVcoeff(file)
load(file,'Vval')
pvar t x_1 x_2 x_3 x_4
Vdeg = 4;
Vmonom = monomials([x_1;x_2;x_3;x_4;t], 0:Vdeg);
[Vcoeff,R,e] = poly2basis(Vval,Vmonom);