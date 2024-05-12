function rok=f1_Element_nodes(de,dn,Nelem,xg,Relem)

Num_grids=size(xg,2);
rok=zeros(Num_grids,1);

xxg=xg(1,:);
yxg=xg(2,:);


for ne=1:Nelem
    
    
    xn1=dn(de(ne,1),1);
    yn1=dn(de(ne,1),2);
    
    xn2=dn(de(ne,2),1);
    yn2=dn(de(ne,2),2);
    
    xn3=dn(de(ne,3),1);
    yn3=dn(de(ne,3),2);
    

    xv=[xn1,xn2,xn3,xn1]; yv=[yn1,yn2,yn3,yn1];
    
    in =inpolygon(xxg,yxg,xv,yv);
    
    
     
    for i=1:Num_grids
        
        if in(i)==1
                                                           
            rok(i)=Relem(ne);

       end
    end
    
    
end

return