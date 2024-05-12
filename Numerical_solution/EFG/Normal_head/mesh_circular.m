close all; clear all;clc

datain_HLD_fine_mesh
% kk=2;
% xe(1,kk:kk+size(nn,2)-1)=nn(1,:);
% xe(2,kk:kk+size(nn,2)-1)=nn(2,:);
% kk=kk+size(nn,2);
% n=0;
% step=5.3;
% for r=step:step:250-step
% teta=0:pi/(3+n):2*pi;
% xe(1,kk:kk+size(teta,2)-1)=r*cos(teta);
% xe(2,kk:kk+size(teta,2)-1)=r*sin(teta);
% kk=kk+size(teta,2);
% n=n+3;
% end


kk=2;
n=0;
step=5;
cc=0;

for r=step:step:250
teta=0:pi/(3+n):2*pi;
xe(1,kk:kk+size(teta,2)-1)=r*cos(teta);
xe(2,kk:kk+size(teta,2)-1)=r*sin(teta);
kk=kk+size(teta,2);
n=n+2.5;
% if r==250
%     cc=cc+1;
%         nn(1,cc)=r*cos(teta);
%         nn(2,cc)=r*sin(teta);
% end
end

x2=xe;


M=size(x2,2);
N=16*4;

x=zeros(2,M+N);

xe(:,1:M)=x2;
xe(:,M+1:M+N)=n_recivNode;




numnod=size(xe,2)
hold on
plot(xe(1,:),xe(2,:),'.k')
axis equal

% cc=0;
% 
% for i=1:numnod
%     r(i)=sqrt((xe(1,i)^2)+(xe(2,i)^2));
%     if (abs(r(i)-250)<0.001)
%         cc=cc+1;
%         nn(:,cc)=xe(:,i);
%     end
% end
% 
% figure
% plot(nn(1,:),nn(2,:),'.')

