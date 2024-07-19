function u=energy_sa(x,h,w,M,c,w_proj,w_hom)



%%%%%%%%% hom_term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 I=rowreshape(x,h,w);
 
 a=diff(I,[],1); % horizontal edges
 b=diff(I,[],2); % vertical edges

 a=a.^2;
 b=b.^2;
 
 border1=(I(:,1)-0); border1=border1.^2;   % PROBE
 border2=(I(:,end)-0); border2=border2.^2; % PROBE
 border3=(I(1,:)-0); border3=border3.^2;   % PROBE
 border4=(I(end,:)-0); border4=border4.^2; % PROBE
 
 %g=zeros(size(I));
 %g(1:end-1,:)=g(1:end-1,:)+a; % horizontal edges
 %g(:,1:end-1)=g(:,1:end-1)+b; % vertical edges
 %hom_term=sum(g(:));
 
 hom_term=sum(a(:))+sum(b(:))+sum(border1)+sum(border2)+sum(border3)+sum(border4); % PROBE
 
%%%%%%%%% hom_term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   

u=0.5*w_proj*(M*x-c)'*(M*x-c)+w_hom*0.5*hom_term;
