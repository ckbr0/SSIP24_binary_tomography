function x=clamp(x,a,b)
%function x=clamp(x,a,b)
% Restrict range of x to [a,b]

x(x<a)=a;
x(x>b)=b;
