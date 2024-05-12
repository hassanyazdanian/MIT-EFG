%%this program make data for FEM code 
%%this program save data(mesh information) in data mfile
              %%%%%% NOTE %%%%%%%%%%%
%%before run this program must export mesh from COMSOL software
%%then save it to program`s folder 

clc;clear;close all

dat=importdata('HN_fine_mesh.mphtxt');
Nnode=str2double(dat.textdata(15,1)); % number of nodes
dn=zeros(Nnode,2);                   % data for nodes
VarName=str2double(dat.textdata(:,1));
Created=dat.textdata(:,2);
by=dat.textdata(:,3);
comsol=dat.textdata(:,4);

R=250;  %the raduis for boundary circle    
for i=1:Nnode
    dn(i,1)=str2double(dat.textdata(i+17,1));   
    dn(i,2)=str2double(dat.textdata(i+17,2));
%     if dn(i,1)^2+dn(i,2)^2> (R-0.001)^2
%         dn(i,3)=-1;
%     end
end

nn=[];
nd=[];



for i=1:size(dn,1)
    r(i)=sqrt(dn(i,1)^2+dn(i,2)^2);
end
for i=1:size(dn,1)
    if (abs(r(i)-R)<=10^-8 )

        i;
        nn=[nn;dn(i,:)];
        nd=[nd i];
    end
end



dn=dn(:,1:2);
%% --------------
for i=1:size(VarName)
    if isequal(cell2mat(Created(i)),'tri')
        j=i;
        Nelem=VarName(i+2,1);  %number of elements
    end
end
de=zeros(Nelem,4);  % data for elements
for i=1:Nelem
     de(i,1)=VarName(i+j+3,1)+1;   %1st node
     de(i,2)=str2double(Created(i+j+3,1))+1; %2nd node
     de(i,3)=str2double(by(i+j+3,1))+1; %3rd node
     de(i,4)=VarName(Nelem+i+j+5,1);
end


% Example=1;  % normal Head
Example=2;   % small deep hemorrhage (HSD)

if Example==2  %  coil design example
  NumSource=32;    
  NumMat=38;
  scale=0.0010;   %the comsol output is in millimeter        
  somat=zeros(NumMat,1);
   
end    

fid = fopen('HN_fine_mesh_data_F1.m','w');
fprintf(fid,'%20s \n ','%		            	 I N dn U de   F I L E' );
fprintf(fid,'%6s \n ','% Date: ');
fprintf(fid,'%15s \n ','% File Number :  ');
fprintf(fid,'%10s \n ','% Subjcet:');
fprintf(fid,'%30s \n ','%.................  << Input Parameters >>   ...................');
fprintf(fid,'%25s \n ','%   nnode  nelem  nsorc  nmat scale ');
fprintf(fid,'%4s %6i %6i  %2i %3i  %3i   %4s \n','dat=[',Nnode,Nelem,NumSource,NumMat,scale,'];');
fprintf(fid,'%30s \n ','%...................... << Nodes >>  ...........................');
fprintf(fid,'%25s \n ','%   Num	   Cons    x      y');
fprintf(fid,'%4s  ',' dn=[');
       
for i=1:Nnode   
   if i~=Nnode
       if i==1
          fprintf(fid,'%12.8f  %12.8f  \n',dn(i,1),dn(i,2));
       else
           fprintf(fid,'%12.8f  %12.8f  \n',dn(i,1),dn(i,2));
       end
   else
      fprintf(fid,'%12.8f  %12.8f %4s \n',dn(i,1),dn(i,2),'];');
   end   

end   

fprintf(fid,'%30s \n',' % ....................... << Elements >> ...........................');
fprintf(fid,'%20s \n',' %   Num.  n1   n2  n3   sorce  mat.');
fprintf(fid,'%4s  ',' de=[');

for i=1:Nelem 
   if i~=Nelem
       if i==1
          fprintf(fid,'%3i  %6i  %6i   %4i \n',de(i,1),de(i,2),de(i,3),de(i,4));
       else
          fprintf(fid,'%10i  %6i  %6i   %4i \n',de(i,1),de(i,2),de(i,3),de(i,4));
       end
   else
     fprintf(fid,'%10i %6i  %6i  %4i %4s \n',de(i,1),de(i,2),de(i,3),de(i,4),'];' );
   end
end 

fprintf(fid,'%30s \n',' % ....................... << Boundary nodes coordinates>> ...........................');
fprintf(fid,'%20s \n',' %   nn.');
fprintf(fid,'%4s  ',' nn=[');

Nvtx=length(nn);

for i=1:Nvtx 
   if i~=Nvtx
       if i==1
          fprintf(fid,'%3i %6i \n',nn(i,1),nn(i,2));
       else
          fprintf(fid,'%10i %6i\n',nn(i,1),nn(i,2));
       end
   else
     fprintf(fid,'%10i %6i %4s\n',nn(i,1),nn(i,2),'];' );
   end
end
fprintf(fid,'%30s \n',' % ....................... <<Boundary nodes>> ...........................');
fprintf(fid,'%20s \n',' %   nd.');
fprintf(fid,'%4s  ',' nd=[');

for i=1:Nvtx 
   if i~=Nvtx
       if i==1
          fprintf(fid,'%3i  \n',nd(1,i));
       else
          fprintf(fid,'%10i\n',nd(1,i));
       end
   else
     fprintf(fid,'%10i %4s\n',nd(1,i),'];' );
   end
end


fprintf(fid,'%30s \n',' % ....................... <<result nodes>> ...........................');
fprintf(fid,'%20s \n',' %   nd.');

 
fprintf(fid,'%20s \n','nn = transpose(nn);');
fprintf(fid,'%20s \n','nd = transpose(nd);');
% fprintf(fid,'%20s \n','RR = transpose(RR);');
%%
fclose(fid);
fprintf('%4s \n',' Terminated.' );
fprintf('%4s \n',' HN_fine_mesh_data_F1.m" ' );
