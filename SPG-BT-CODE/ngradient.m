% NUMERICAL GRADIENT
% Tibor Lukic 2014

% input: 
% FUN(x) - functional
% xpoint - vector

%  calculate numerical gradient of the function FUN(x) at the point xpoint

xpoint=[1 1 1]';

n=length(xpoint);

delta_step=10^(-4);

v=ones(n,1);

X=delta_step*diag(v);

ngradient=zeros(n,1);

for i=1:n,
    
    f1=FUN(xpoint-X(i,:));
    f2=FUN(xpoint+X(i,:));
    ngradient(i)=(f2-f1)/(2*delta_step);
    
end;