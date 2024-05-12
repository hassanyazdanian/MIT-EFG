% function [AFE]=EFG_FE(rho_EFG,rho_FE,NC,current,Electrodes_segments,Nnode,Nelem,dn,de,db,gp,Base,numnod,Ns,ninter,xe,dm,b , c , DELTA)
close all;clear all;clc
HN_fine_mesh_data_F1

data_FE1
ndivl1 = 180;ndivw1 =180;
ndivl = 180;ndivw =180;
dmax=1.5;    %scaling parameter between 2-4
ndivlq =65;ndivwq =65;  % SET UP QUADRATURE CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    FE part   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nnode=length(dn(:,1));

Nelem=length(de(:,1));

nsorc=dat(1,3);% Number of sources
nmat=dat(1,4);% Number of materials
scale=dat(1,5);% Scale
jj=sqrt(-1);   %imaginary symbol

dn=dn*scale;

ncoil=16;
L=0.2;     % Length of coils (Z Direction) 20 cm
D=0;       % Width of coils (5.8) are negligible respect to their lengths
Turn=10;    % Turn in each coil

Acoil=(5*scale*30*scale);
magsorc=ds(:,2)./Acoil;
angsorc=pi/180*ds(:,3);



excit=zeros(ncoil,nmat);
mm=1;

for ic=1:ncoil
    for jc=1:2
        excit(ic,coil(ic,jc))=magsorc(mm)*cos(angsorc(mm))+jj*magsorc(mm)*sin(angsorc(mm));
        mm=mm+1;
    end
end


% recivNode=zeros(ncoil,4);
% 
% for i=1:ncoil
%     for j=1:4
%         for k=1:size(dn,1)
%             if ((abs(dn(k,1)-x_recivNode(i,j))<eps)&&(abs(dn(k,2)-y_recivNode(i,j))<eps))
%                 recivNode(i,j)=k;
%             end
%         end
%     end
% end
recivNode=[1614,2199,1962,2459;2273,2529,2348,2656;2367,2658,2356,2660;2252,2617,2132,2491;2021,2494,1906,2410;1781,2327,1655,2115;1542,2112,1424,1986;1310,1875,1194,1518;1000,1517,882,1191;601,1088,493,779;320,771,237,480;129,478,198,524;293,564,305,760;461,851,544,1054;732,1164,1143,1499;1362,1727,1592,2191];
% x_recivNode=zeros(size(recivNode));
% y_recivNode=zeros(size(recivNode));
% n_recivNode=zeros(2,64);
% 
% % plot receive nodes
% figure;
% hold on
% coun=0;
% for i1=1:ncoil
%     for j=1:4
%         coun=coun+1;
%         x_recivNode(i1,j)=dn(recivNode(i1,j),1);
%         y_recivNode(i1,j)=dn(recivNode(i1,j),2);
%         n_recivNode(1,coun)=dn(recivNode(i1,j),1);
%         n_recivNode(2,coun)=dn(recivNode(i1,j),2);
%         plot(x_recivNode(i1,j),y_recivNode(i1,j),'or')
%     end
% end

[b , c , Area]=SeMatrix(dn,de);
exper=0;                              %experiment counter

K=sparse(Nnode,Nnode);

for ie=1:Nelem %{10}
    for m=1:3
        for n=1:3
            if m==n
                delta=1;
            else
                delta=0;
            end
            
            K(de(ie,m),de(ie,n))= K(de(ie,m),de(ie,n))+(1/(4*Area(ie,1)))*(1/nu0)*((b(ie,m).*b(ie,n)+c(ie,n).*c(ie,m)));
%             +(Area(ie,1)/12)*(1i*w*dmm1(de(ie,4),3))*(1+delta);
            
        end
    end
end %{10}

   

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    EFG part   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% DEFINE BOUNDARIES/PARAMETERS
Lx = 500;    %length
Ly = 500;     %height 


xspac = Lx/ndivl;
yspac = Ly/ndivw;
Nc=length(nc);
% nc(1:2,:)=nc(1:2,:)*scale;
[xe,conn1,conn2] = mesh1(Lx,Ly,ndivl1,ndivw1,nc);
xe=xe*scale;

% NODE=[dn;xe'];
% pdeplot(dn',[],de')
% hold on
% plot(xe(1,:),xe(2,:),'.b')
numnod=size(xe,2);
% Nnode+numnod-Nc

dm = sparse(2,numnod);     %The size of the domain of influence at a node
dm(1,1:numnod)=dmax*xspac*scale*ones(1,numnod);
dm(2,1:numnod)=dmax*yspac*scale*ones(1,numnod);

[xc,conn,numcell,numq] = mesh2(Lx,Ly,ndivlq,ndivwq);


% SET UP GAUSS POINTS, WEIGHTS, AND JACOBIAN FOR EACH CELL
quado = 4;
[gauss] = pgauss(quado);
numq2 = numcell*quado^2;
gp = sparse(4,numq2);

[gp] = egauss(xc,conn,gauss,numcell,nc);
gp=gp*scale;

Nelem1=length(de1(:,1));

for ie=1:Nelem1 %{10}
    sige(ie)=dmm1(de1(ie,4),3);   
end %{10}
dn1=dn1*scale;
Sig=f1_Element_nodes(de1,dn1,Nelem1,gp,sige);


%  sign=pdeprtni(dn1',de1',sige);
%  Sig = griddata(dn1(:,1),dn1(:,2),sign,gp(1,:),gp(2,:),'cubic');


% figure
% plot3(gp(1,:),gp(2,:),Sig,'.')
% 
% fghfghjk


KEFG=sparse(numnod,numnod);
BEFG = sparse(numnod,1);

% LOOP OVER GAUSS POINTS TO ASSEMBLE DISCRETE EQUATIONS
 for j=1:size(gp,2)
   gg=gp(:,j);

   gpos = gg(1:2);
   weight = gg(3);
   jac = gg(4);


% DETERMINE NODES IN NEIGHBORHOOD OF GAUSS POINT
   en = domain(gpos,xe,dm,numnod);
   LL = length(en);
   [phi,dphix,dphiy] = shape(gpos,dmax,xe,en,dm);

   KEFG(en,en) = KEFG(en,en)+(weight*(1/nu0)*jac).*(dphix'*dphix+dphiy'*dphiy)+(weight*(1i*w*Sig(j))*jac).*(phi'*phi);

end




ninter=nc(:,1:end-1);
ninter(1:2,:)=ninter(1:2,:)*scale;
Ns=length(ninter);
GFE = sparse(Nnode,Ns);

for i=1:(Ns)
   
  GFE(ninter(3,i),i)=GFE(ninter(3,i),i)+1;  %Eq32
            
end
GEFG = sparse(numnod,Ns);
for i=1:(Ns)
   
  gpos = ninter(1:2,i);
  
   v = domain(gpos,xe,dm,numnod);

    [phi,dphix,dphiy] = shape(gpos,dmax,xe,v,dm);
   LL = length(v);
   for n=1:LL
      GEFG(v(n),i)=GEFG(v(n),i)-phi(n);%Eq33
   end
           
end
for ic=1:ncoil
% for ic=1
    A=sparse(Nnode,1);
    B=sparse(Nnode,1);
    fe(:,1)=sparse(excit(ic,de(:,4)));
    for ie=1:Nelem
        %f is a excitation source
        for m=1:3
            Be(ie,m)=Area(ie)*fe(ie)/3;            %Be is a element excitation matrix
            B(de(ie,m),1)=B(de(ie,m),1)+Be(ie,m);
        end
    end
    KFE=(K);
    BFE=(B);
    % Knew(nd,nd)=10^70.*eye(size(nd,2));
    % Bnew(nd)=10^70.*AT;
    
    Q=sparse(Nnode,1);  %Dirichlet Boundary Condition
    
    KFE(nd,:)=sparse(size(nd,2),size(KFE,2));
    KFE(:,nd)=sparse(size(KFE,2),size(nd,2));
    KFE(nd,nd)=eye(size(nd,2));
    Q(nd,1)=0;
    BFE=B-K*Q;
    BFE(nd)=Q(nd);
    m=([KFE zeros(Nnode,numnod)  GFE;zeros(numnod,Nnode) KEFG  GEFG;GFE.' GEFG.' zeros(Ns,Ns)]);
    f=[BFE;BEFG;zeros(Ns,1)];
    A1=m\f;

    AFE=A1(1:Nnode);
%     Aefg=A1(Nnode+1:Nnode+numnod);
%     
%     parfor ind=1:length(xe)
%     gg=xe(:,ind);
%     gpos = gg(1:2);
%     v = domain(gpos,xe,dm,numnod);
%     [phi,dphix,dphiy] = shape(gpos,dmax,xe,v,dm);
%      AEFG(ind,1) = phi*Aefg(v,1);
%     end
%     A=full([AFE;AEFG]);
    
    
  
    Aind(ic,:)=AFE.';%induction voltage vector for each excit
    for ir=ic+1:ncoil
        exper=exper+1;
        V_F1_HLD(exper,1)=-jj*w*Turn*L*(mean([Aind(ic,(recivNode(ir,1))),Aind(ic,(recivNode(ir,2)))])-mean([Aind(ic,(recivNode(ir,3))),Aind(ic,(recivNode(ir,4)))]));%measuring voltages
        % warning: the sign must be minus!
        % it must be noted that the minuse sign in voltage calculation it is
        % for correction of receiving nodes direction. In fact receiving node 2
        % must be minus from receiving node 1.
    end
end
V_F1_HLD_FE=[-0.209384556572310 + 119.363944026187i;0.0136058214423580 + 34.6003791625659i;0.0751261138466377 + 14.2792662240434i;0.0885013637830983 + 6.86269631522702i;0.0899946520125323 + 3.70729287114460i;0.0921069847846824 + 2.27097747726567i;0.0961896977853103 + 1.63651178407778i;0.0976549255995654 + 1.45394353602588i;0.0939314050413064 + 1.63405338924460i;0.0860115059755323 + 2.27090092891099i;0.0771566521790619 + 3.70950284861134i;0.0680352744743281 + 6.87561720753873i;0.0471492596516829 + 14.3021595755245i;-0.0234788500847273 + 34.6678435315867i;-0.236882942135318 + 119.810917569208i;-0.295158753541598 + 119.819375102486i;-0.124813204117370 + 34.6496802168630i;-0.0310203444789236 + 14.3116056590258i;0.0188221546568140 + 6.87659287418160i;0.0507445718645771 + 3.71030061237353i;0.0750067793679523 + 2.27172297583312i;0.0932168610059753 + 1.63479298487878i;0.104081461609274 + 1.45213441814089i;0.109229560042086 + 1.63414147316137i;0.117151946671274 + 2.26853380179440i;0.134005248737534 + 3.70471183027095i;0.158073394512076 + 6.85564498068537i;0.169608370715176 + 14.2785814935798i;0.0822888424230955 + 34.6538500941470i;-0.282584582009360 + 119.718517015473i;-0.150791711077644 + 34.6364444955389i;-0.0584855097551441 + 14.3315550600829i;0.00550438665831167 + 6.88605967097467i;0.0522796801605873 + 3.71297907159453i;0.0864735671381426 + 2.27124801806457i;0.106741022787758 + 1.63410275535608i;0.114567961516523 + 1.45330509615359i;0.122992340230197 + 1.63405655104552i;0.140742682466496 + 2.26877245634865i;0.169944515841236 + 3.69779205387430i;0.201652336287182 + 6.85004484160055i;0.182516365908367 + 14.2841744022361i;-0.241143980878981 + 119.629058448122i;-0.141868761428354 + 34.6060807105590i;-0.0488775782107427 + 14.3344007829154i;0.0265567394820826 + 6.88150102610811i;0.0818843073062670 + 3.70927620746938i;0.112258125031102 + 2.26897812850893i;0.118792704114273 + 1.63476189740753i;0.120784487054302 + 1.45341635598902i;0.129374341929172 + 1.63551942455893i;0.148059398932823 + 2.26590839347427i;0.172386215779774 + 3.69689648802076i;0.169765057725555 + 6.85813707540055i;-0.233483353435458 + 119.602880041678i;-0.130503429099223 + 34.6355182265933i;-0.0162135691035225 + 14.3237456525833i;0.0750099163019823 + 6.87254219816319i;0.123925073890069 + 3.70365158729438i;0.130681364735666 + 2.26829118930226i;0.125230416026797 + 1.63419135756073i;0.123265927848168 + 1.45508083167333i;0.129909391197086 + 1.63397585021261i;0.143276861780518 + 2.26632114110496i;0.143129255855062 + 3.70425539159847i;-0.247640716373398 + 119.632552613212i;-0.102013683796664 + 34.6259557728859i;0.0495340891029338 + 14.3026418162153i;0.135608629880285 + 6.85860913451990i;0.149060153079386 + 3.70006757631967i;0.137623927270481 + 2.26647266554405i;0.125596640399933 + 1.63598425688263i;0.121103696632171 + 1.45367742961989i;0.124368821728475 + 1.63351179292415i;0.122757524119003 + 2.26928798369980i;-0.262582353302203 + 119.686676156012i;-0.0353662609353900 + 34.6930749013405i;0.123736145405234 + 14.2990686181944i;0.161670051555655 + 6.85871710628239i;0.151244104888675 + 3.69936048646121i;0.132332331206411 + 2.26959017474199i;0.119106994077537 + 1.63498437026403i;0.114573493619646 + 1.45322319887019i;0.111404000196198 + 1.63569390850811i;-0.233898682129891 + 119.853681744142i;0.0293449096245762 + 34.6728607787403i;0.132134597837236 + 14.3069599167729i;0.144515424072066 + 6.86649778514138i;0.129845061244167 + 3.70980209666574i;0.114493626787013 + 2.27252661256844i;0.106820610726655 + 1.63646310043918i;0.104601863668510 + 1.45581656875717i;-0.224363591758917 + 119.245479281542i;-0.0130903518799523 + 34.6366039351267i;0.0704595967605629 + 14.2993752801056i;0.0886961900828049 + 6.86990284178731i;0.0880022778704640 + 3.70768382857501i;0.0876854839419996 + 2.26949087764138i;0.0930732205227278 + 1.63487161340007i;-0.258218306932490 + 119.557182742365i;-0.0822258740678146 + 34.6793087538144i;-0.000397992499887689 + 14.3200634522721i;0.0335544025252653 + 6.87766957686541i;0.0527480429035873 + 3.70784678422263i;0.0733546153938874 + 2.26881079562997i;-0.245914351690408 + 119.586689009086i;-0.121110196800770 + 34.6330111699760i;-0.0452018473145437 + 14.3384502893157i;0.00238245724957668 + 6.88470033247971i;0.0455391335724354 + 3.70995095840530i;-0.233502330037848 + 119.371745082548i;-0.143603234158837 + 34.6506605266039i;-0.0674354364631769 + 14.3437813106121i;0.00743240846178545 + 6.88060726183219i;-0.238816669927077 + 118.875341707735i;-0.159213536809978 + 34.6292042617008i;-0.0467743309290061 + 14.3187723421736i;-0.274222655271636 + 119.093726324379i;-0.136716537320368 + 34.6281394871194i;-0.289139124916458 + 119.403076285452i];


plot(real(V_F1_HLD),'r')
hold on
plot(real(V_F1_HLD_FE),'b')

figure

plot(imag(V_F1_HLD),'r')
hold on
plot(imag(V_F1_HLD_FE),'b')

AFE=full(AFE);
% figure
% pdeplot(dn',[],de','xydata',real(AFE))
% figure
% pdeplot(dn',[],de','xydata',imag(AFE))
