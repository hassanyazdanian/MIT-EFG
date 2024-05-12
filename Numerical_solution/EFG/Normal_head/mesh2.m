%-------------------------------------------------------------------------
% MESH GENERATION PROGRAM FOR 2D BEAM IN BENDING
function[x,conn,numcell,numq] = mesh2(length,height,ndivl,ndivw)

% INPUT DATA
numcell= ndivw*ndivl;
numq = (ndivl+1)*(ndivw+1); 

% SET UP NODAL COORDINATES
for i = 1:(ndivl+1)
   for j = 1:(ndivw+1)
   x(1,((ndivw+1)*(i-1) +j))= (length/ndivl)*(i-1)-length/2;
   x(2,((ndivw+1)*(i-1) +j))= -(height/ndivw)*(j-1)+height/2;
   end
end

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
