close all;clear all;clc

syms a [14 1]
syms sigma [1 8]
syms n
syms phi
syms r [1 7]
syms alpha   %alpha=j*w*mu
syms J
syms mu
syms eq [14 1]
syms S [14 1]
syms A [52 1]
syms I
syms w

for i=1:7
    if i==1
        eqs(i,:)=a(i)*A(1)-a(i+1)*A(2)-a(i+2)*A(3)==0;
        eqs(i+7,:)=a(i)*A(27)-a(i+1)*A(28)-a(i+2)*A(29)==0;
        m=0;
        nn=3;
    elseif (i~=7)
        eqs(i,:)=a(i+m)*A(nn+1)+ a(i+m+1)*A(nn+2)-a(i+m+2)*A(nn+3)-a(i+m+3)*A(nn+4)==0;
        eqs(i+7,:)=a(i+m)*A(nn+27)+ a(i+m+1)*A(nn+28)-a(i+m+2)*A(nn+29)-a(i+m+3)*A(nn+30)==0;
        m=m+1;
        nn=nn+4;
    elseif i==7
        eqs(i,:)=a(i+m)*A(nn+1)+ a(i+m+1)*A(nn+2)-a(i+m+2)*A(nn+3)==0;
        eqs(i+7,:)=a(i+m)*A(nn+27)+ a(i+m+1)*A(nn+28)-a(i+m+2)*A(nn+29)==I;
    end
    
end

            
             
    
% [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14]=solve(eqs,a);

% [A,b] = equationsToMatrix(eqs,a);
% z = A\b

% simplify(eq)

