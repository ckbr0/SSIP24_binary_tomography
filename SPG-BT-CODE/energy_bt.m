function u=energy_bt(x,h,w,M,c,w_hom,mu)


%%%%%%%%% hom_term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 I=rowreshape(x,h,w);
 a=diff(I,[],1); % horizontal edges
 b=diff(I,[],2); % vertical edges

 a=a.^2;
 b=b.^2;
 g=zeros(size(I));
 g(1:end-1,:)=g(1:end-1,:)+a; % horizontal edges
 g(:,1:end-1)=g(:,1:end-1)+b; % vertical edges
 hom_term=sum(g(:));
 
%%%%%%%%% hom_term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=ones(h*w,1);  

u=0.5*(M*x-c)'*(M*x-c)+w_hom*0.5*hom_term+0.5*mu*(x'*(E-x));

% u=0.5*(M*x-c)'*(M*x-c)+w_hom*0.5*hom_term+0.5*mu*(x'*(E-x))^2;  % probe
