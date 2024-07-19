function ngradient=num_gradient_par(FUN,xpoint)

% calculate numerical gradient of the function FUN at the point xpoint
% Tibor Lukic 2015

% input: 
% FUN - function/functional, string type
% xpoint - vector

n=length(xpoint);
delta_step=10^(-3); % suggested value: -3

% X=delta_step*diag(ones(n,1));
ngradient=zeros(n,1);
x0=ngradient;


parfor i=1:n,   % probe: parfor
    
    
    xx=x0;
    xx(i)=delta_step;
    
      f1=orientation_vec(xpoint-xx); % PROBE
      f2=orientation_vec(xpoint+xx); % PROBE
    
     % f1=feval( FUN, xpoint-xx); 
     % f2=feval( FUN, xpoint+xx); 
    
    ngradient(i)=(f2-f1)/(2*delta_step);   
    
end;