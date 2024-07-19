function ht=hom_term(I)
 
 a=diff(I,[],1); % horizontal edges
 b=diff(I,[],2); % vertical edges

 a=a.^2;
 b=b.^2;
 
 border1=(I(:,1)-0); border1=border1.^2;   % PROBE
 border2=(I(:,end)-0); border2=border2.^2; % PROBE
 border3=(I(1,:)-0); border3=border3.^2;   % PROBE
 border4=(I(end,:)-0); border4=border4.^2; % PROBE
 

 
ht=sum(a(:))+sum(b(:))+sum(border1)+sum(border2)+sum(border3)+sum(border4); % PROBE