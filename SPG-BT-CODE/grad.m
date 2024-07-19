% Image gradient calculation
function p=grad(I)

  a=diff(I,[],1); % horizontal edges
  b=diff(I,[],2); % vertical edges


  a=a.^2;
  b=b.^2;
  % gradient
  g=zeros(size(I));
  g(1:end-1,:)=g(1:end-1,:)+a; % horizontal edges
  g(:,1:end-1)=g(:,1:end-1)+b; % vertical edges
  

  p=sum(g(:));