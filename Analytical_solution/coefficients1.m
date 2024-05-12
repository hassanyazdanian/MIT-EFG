close all;clear all;clc

syms sigma [1 8]
syms n
syms phi
syms r [1 7]
syms alpha   %alpha=j*w*mu
syms J
syms mu
syms S [14 1]
syms A [14 14]
syms I [14 1]

I(14)=mu*J;


for i=1:7
    if i==1
        A(1,1)=besseli(n,r(i)*(alpha*sigma(i))^(1/2));
        A(1,2)=-besseli(n, r(i)*(alpha*sigma(i+1))^(1/2));
        A(1,3)=-besselk(n,r(i)*(alpha*sigma(i+1))^(1/2));
        A(8,1)=(alpha*sigma(i))^(1/2)*(besseli(n-1,r(i)*(alpha*sigma(i))^(1/2)) + besseli(n+1,r(i)*(alpha*sigma(i))^(1/2)))/2;
        A(8,2)=-((alpha*sigma(i+1))^(1/2)*(besseli(n-1,r(i)*(alpha*sigma(i+1))^(1/2))+besseli(n+1,r(i)*(alpha*sigma(i+1))^(1/2)))/2);
        A(8,3)=(alpha*sigma(i+1))^(1/2)*(besselk(n-1,r(i)*(alpha*sigma(i+1))^(1/2))+besselk(n + 1,r(i)*(alpha*sigma(i+1))^(1/2)))/2;
        nn=0;
    elseif (i~=6)&&(i~=7)
        A(i,nn+i)=besseli(n,r(i)*(alpha*sigma(i))^(1/2));
        A(i,nn+i+1)=besselk(n,r(i)*(alpha*sigma(i))^(1/2));
        A(i,nn+i+2)=-besseli(n,r(i)*(alpha*sigma(i+1))^(1/2));
        A(i,nn+i+3)=-besselk(n, r(i)*(alpha*sigma(i+1))^(1/2));
        A(i+7,nn+i)=(alpha*sigma(i))^(1/2)*(besseli(n-1,r(i)*(alpha*sigma(i))^(1/2))+besseli(n+1,r(i)*(alpha*sigma(i))^(1/2)))/2;
        A(i+7,nn+i+1)=-(alpha*sigma(i))^(1/2)*(besselk(n-1,r(i)*(alpha*sigma(i))^(1/2))+besselk(n+1,r(i)*(alpha*sigma(i))^(1/2)))/2;
        A(i+7,nn+i+2)=-((alpha*sigma(i+1))^(1/2)*(besseli(n-1,r(i)*(alpha*sigma(i+1))^(1/2))+besseli(n+1,r(i)*(alpha*sigma(i+1))^(1/2)))/2);
        A(i+7,nn+i+3)=+((alpha*sigma(i+1))^(1/2)*(besselk(n-1,r(i)*(alpha*sigma(i+1))^(1/2))+besselk(n+1,r(i)*(alpha*sigma(i+1))^(1/2)))/2);
        nn=nn+1;
    elseif i==6
        A(i,nn+i)=besseli(n,r(i)*(alpha*sigma(i))^(1/2));
        A(i,nn+i+1)=besselk(n,r(i)*(alpha*sigma(i))^(1/2));
        A(i,nn+i+2)=-r(i);
        A(i,nn+i+3)=-r(i)^(-1);
        A(i+7,nn+i)=(alpha*sigma(i))^(1/2)*(besseli(n-1, r(i)*(alpha*sigma(i))^(1/2))+besseli(n+1,r(i)*(alpha*sigma(i))^(1/2)))/2;
        A(i+7,nn+i+1)=-(alpha*sigma(i))^(1/2)*(besselk(n-1,r(i)*(alpha*sigma(i))^(1/2))+besselk(n + 1,r(i)*(alpha*sigma(i))^(1/2)))/2;
        A(i+7,nn+i+2)=-1;
        A(i+7,nn+i+3)=r(i)^(-2);
        nn=nn+1;
        elseif i==7
        A(i,nn+i)=r(i);
        A(i,nn+i+1)=r(i)^(-1);
        A(i,nn+i+2)=-r(i)^(-1);
        A(i+7,nn+i)=1;
        A(i+7,nn+i+1)=-r(i)^(-2);
        A(i+7,nn+i+2)=r(i)^(-2);

    end
    
end


S=A\I;
