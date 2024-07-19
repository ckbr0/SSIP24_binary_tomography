function ngradient=num_gradient(FUN,xpoint) %#codegen

% calculate numerical gradient of the function FUN at the point xpoint
% Tibor Lukic 2015

% input: 
% FUN - function/functional, string type
% xpoint - vector

n=length(xpoint);
delta_step=10^(-3);
f1=0; f2=0;


ngradient=zeros(n,1);
i=1;
while i<(n+1), 
    x_back=xpoint;
    x_back(i)=x_back(i)-delta_step;
    
    x_forward=xpoint;
    x_forward(i)=x_forward(i)+delta_step;
    
    f1=feval( FUN, x_back );
    f2=feval( FUN, x_forward );
    ngradient(i)=(f2-f1)/(2*delta_step);   
    i=i+1;
end;