%--------------------------------------------------------------------------
% EFG shape function and it's derivative with linear base
function  [phi,dphix,dphiy] = shape(gpos,dmax,x,v,dm)

L = length(v);
won = ones(1,L);
nv = x(1:2,v);
p = [won;nv];
dif = gpos*won-nv;
t = dm(1:2,v)/dmax;

% WEIGHTS--W AND DW ARE VECTORS
[w,dwdx,dwdy] = cubwgt(dif,t,v,dmax,dm);
B = p.*[w;w;w];
pp = zeros(3);
aa = zeros(3);
daax = zeros(3);
daay = zeros(3);
for i=1:L
   pp = p(1:3,i)*p(1:3,i)';
   aa = aa+w(1,i)*pp;
   daax = daax+dwdx(1,i)*pp;
   daay = daay+dwdy(1,i)*pp;
end
pg = [1 gpos'];
[L,U,PERM] = lu(aa);
for i=1:3
   if i==1
      C = PERM*pg';
   elseif i==2
      C = PERM*([0 1 0]' - daax*gam(1:3,1));
   elseif i==3
      C = PERM*([0 0 1]' - daay*gam(1:3,1));
   end
D1 = C(1);
D2 = (C(2) - L(2,1)*D1);
D3 = (C(3) - L(3,1)*D1 - L(3,2)*D2);
gam(3,i) = D3/U(3,3);
gam(2,i) = (D2 - U(2,3)*gam(3,i))/(U(2,2));
gam(1,i) = (D1 - U(1,2)*gam(2,i)-U(1,3)*gam(3,i))/U(1,1);
end

phi = gam(1:3,1)'*B;
dbx = p.*[dwdx;dwdx;dwdx];
dby = p.*[dwdy;dwdy;dwdy];
dphix = gam(1:3,2)'*B + gam(1:3,1)'*dbx;
dphiy = gam(1:3,3)'*B + gam(1:3,1)'*dby;

