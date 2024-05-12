%-------------------------------------------------------------------------
function v = domain(gpos,x,dm,nnodes)
%DETERMINES NODES IN DOMAIN OF INFLUENCE OF POINT GPOS
dif = gpos*ones(1,nnodes)-x;
a = dm-abs(dif);
b = [a;zeros(1,nnodes)];
c = all(b>=-100*eps);
v = find(c);
