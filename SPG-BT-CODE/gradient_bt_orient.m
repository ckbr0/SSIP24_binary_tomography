function u=gradient_bt_orient(x,h,w,M,c,w_orient,mu)

global orient_given;

E=ones(h*w,1); 


I=rowreshape(x,h,w); orient_x=orientation(I);

if isnan(orient_x),
    g=gradient_grad_vector(x,h,w);
else
    g=(orient_x-orient_given)*gradient_grad_vector(x,h,w);
end;
    

u=M'*(M*x-c)+w_orient*g-mu*(x-0.5*E);