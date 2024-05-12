%-------------------------------------------------------------------------
% MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING
function[x,conn,numcell] = mesh1(length,height,ndivl,ndivw,nc)

nn(1,:)=nc(1,:);
nn(2,:)=nc(2,:);
nn(3,:)=nc(3,:);
Ns=size(nn,2);
numcell= ndivw*ndivl;
kk=0;
% SET UP NODAL COORDINATES
for i = 1:(ndivl+1)
   for j = 1:(ndivw+1)
   x1(1,((ndivw+1)*(i-1) +j))= (length/ndivl)*(i-1)-length/2;
   x1(2,((ndivw+1)*(i-1) +j))= -(height/ndivw)*(j-1)+height/2;

   end
end

kk=0;

for i=1:size(x1,2)
    
    in =inpolygon(x1(1,i),x1(2,i),nc(1,:),nc(2,:));
    
    if (in==1)
        kk=kk+1;
        x2(:,kk)= x1(:,i);
    end
end


M=size(x2,2);
N=size(nn,2);
% O=size(xg,2);
x=zeros(2,M+N);

x(:,1:M)=x2;
x(:,M+1:M+N)=nn(1:2,:);
% for i=1:N
%     x(:,M+i)=nn(1:2,i);
% end
%  
% for i=1:O
%     x(:,M+N+i)=xg(1:2,i);
% end
 O=0;
TF=0;

while TF==0
  for i=1:M+N+O-1
    if (x(1,i)>x(1,i+1))
        k1=x(1,i);
        k2=x(2,i);
        x(:,i)=x(:,i+1);
        x(1,i+1)=k1;
        x(2,i+1)=k2;
    end
  end
  TF = issorted(x(1,:));
end
  for i=1:M+N+O-1
          if (abs(x(1,i)-x(1,i+1))<10^-010)  
             if (x(2,i)<x(2,i+1))
                  k1=x(2,i);
                  x(2,i)=x(2,i+1); 
                  x(2,i+1)=k1;
             end
          end
  end
       

    
numq=size(x,2);

% SET UP CONNECTIVITY ARRAY
for j = 1:ndivl
   for i = 1:ndivw
   elemn = (j-1)*ndivw + i;
   nodet(elemn,1) = elemn + (j-1);
   nodet(elemn,2) = nodet(elemn,1) + 1;
   nodet(elemn,3) = nodet(elemn,2)+ndivw+1;
   nodet(elemn,4) = nodet(elemn,3)-1;
   end
end
conn = nodet';
