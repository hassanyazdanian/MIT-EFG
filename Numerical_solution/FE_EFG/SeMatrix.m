function  [b , c , DELTA]=SeMatrix(P,T)


Nelem=length(T(:,1));

%%save cordinate of nodes for each element at x & y matrix
x=zeros(Nelem,3); 
y=zeros(Nelem,3);
for i=1:Nelem
    for j=1:3
   x(i,j)= P(T(i,j),1);  %coordinate of nodes(x)
   y(i,j)= P(T(i,j),2);  %coordinate of nodes(y)
    end 
end
%%%Calculation of parameters of nodes for each element
c=zeros(Nelem,3);
b=zeros(Nelem,3);   


for i=1:Nelem

        b(i,1)=y(i,2)-y(i,3);
        c(i,1)=x(i,3)-x(i,2);
        

        b(i,2)=y(i,3)-y(i,1);
        c(i,2)=x(i,1)-x(i,3);
        

        b(i,3)=y(i,1)-y(i,2);
        c(i,3)=x(i,2)-x(i,1);
end
DELTA=zeros(Nelem,1);    
for i=1:Nelem
    DELTA(i,1)=1/2*(b(i,1).*c(i,2)-b(i,2).*c(i,1));
end
