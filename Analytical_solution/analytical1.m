close all;clear ;clc

sigma=[0.8624 0.1582 0.2917 2.002 0.0828 0.6168 0 0];
r=0.001*[30 54 62 66 74 80 100];
mu=4*pi*10^-7;
w=2*pi*1e7;
alpha=1i*w*mu;
J=1;
N=1;
% I=mu*J*sin(phi_source);
R_step=0.1;
phi_step=pi/2;
u=zeros(1,N);

I=mu*J;
for n=1:N
    [S1(n),S2(n),S3(n),S4(n),S5(n),S6(n),S7(n),S8(n),S9(n),S10(n),S11(n),S12(n),S13(n),S14(n)]=f1(sigma,r,w,alpha,I,n);
end

indx=0;  % index


for R=0:R_step:r(end)    % varying from Head center to the Scalp
    %for phi=0:phi_step:2*pi-phi_step
    for phi=pi/2
        for n=1:N
            if (R>=0 && R<=r(1))
                u(n)=S1(n)*besseli(n,R*(alpha*sigma(1))^(1/2))*sin(n*phi); %u1
            elseif (R>r(1) && R<=r(2))
                u(n)=(S2(n)*besseli(n,R*(alpha*sigma(2))^(1/2))+S3(n)*besselk(n,R*(alpha*sigma(2))^(1/2)))*sin(n*phi);  %u2
            elseif (R>r(2) && R<=r(3))
                u(n)=(S4(n)*besseli(n,R*(alpha*sigma(3))^(1/2))+S5(n)*besselk(n,R*(alpha*sigma(3))^(1/2)))*sin(n*phi);  %u3
            elseif (R>r(3) && R<=r(4))
                u(n)=(S6(n)*besseli(n,R*(alpha*sigma(4))^(1/2))+S7(n)*besselk(n,R*(alpha*sigma(4))^(1/2)))*sin(n*phi);  %u4
            elseif (R>r(4) && R<=r(5))
                u(n)=(S8(n)*besseli(n,R*(alpha*sigma(5))^(1/2))+S9(n)*besselk(n,R*(alpha*sigma(5))^(1/2)))*sin(n*phi);  %u5
            elseif (R>r(5) && R<=r(6))
                u(n)=(S10(n)*besseli(n,R*(alpha*sigma(6))^(1/2))+S11(n)*besselk(n,R*(alpha*sigma(6))^(1/2)))*sin(n*phi);  %u6
            elseif (R>r(6) && R<=r(7))
                u(n)=(S12(n)*R+S13(n)*R^(-1))*sin(phi);  %u7
            elseif R>r(7)
                u(n)=S14(n)*R^(-1)*sin(phi);  %u8
            else
                error('Error: please insert correct radii')
            end
        end
    end
    indx=indx+1;
    uT(indx)=sum(u);
    u=[];
end
plot(abs(uT))
