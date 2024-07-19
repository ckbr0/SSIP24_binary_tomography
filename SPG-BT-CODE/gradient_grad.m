% Gradient of the gradient term

function g=gradient_grad(I)

% I image in matrix form


  
a=diff(I,[],1); % vertical edges : x_{i+w}-x_i
b=diff(I,[],2); % horizontal edges : x_{i+1}-x_i
  
g=zeros(size(I));
g(2:end,:)=g(2:end,:)+a;     %edges above
g(1:end-1,:)=g(1:end-1,:)-a; %edges below

g(:,2:end)=g(:,2:end,:)+b;   %edges to the left
g(:,1:end-1)=g(:,1:end-1)-b; %edges to the right

g(:,1)=g(:,1)-I(:,1);         % PROBE
g(:,end)=g(:,end)-I(:,end);   % PROBE
g(1,:)=g(1,:)-I(1,:);         % PROBE
g(end,:)=g(end,:)-I(end,:);   % PROBE

