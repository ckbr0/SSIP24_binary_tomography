function u=fun_probe(x)

% calculate numerical gradient of the function FUN at the point xpoint
% Tibor Lukic 2014

% input: 
% FUN - function/functional, string type
% xpoint - vector

n=length(xpoint);
delta_step=10^(-3);
X=delta_step*diag(ones(n,1));
ngradient=zeros(n,1);
i=1;

while i<(n+1),  
    f1=feval( FUN, xpoint-X(:,i) );
    f2=feval( FUN, xpoint+X(:,i) );
    ngradient(i)=(f2-f1)/(2*delta_step);   
    i=i+1;
end;


