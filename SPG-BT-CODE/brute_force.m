function out_im=brute_force(A,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRUT FORCE algorithm
% for Binary Tomography  on Hexagonag grid
% Ref. Herman and Kuba: Discrete Tomography, Springer, 1999. 

% A - projection matrix
% b - projection vector
% out_im - free pixels are marked by 0.5 intensity 
% Tibor Lukic, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(sprintf('BRUTE FORCE ALGORITHM'));

 % A=sparse(A); b=sparse(b);
 
 n_rays=size(A,1);
 n_pixels=size(A,2);
 x_current=zeros(n_pixels,1); % Every pixel is black
 free_pixels=logical( x_current==x_current ); % All pixels are free 
 i_count=1;
 
while i_count<100,

 for i=1:n_rays,
     
     % STEP ONE
     
     ray_row=A(i,:);
     ray_pos=logical(ray_row~=0);
        
     if ray_row*x_current==b(i),
                  
         x_zeros=zeros(size(free_pixels));
         x_zeros(ray_pos)=free_pixels(ray_pos);
         x_zeros=logical(x_zeros)'; % free pixels on the line
         free_pixels(x_zeros)=0; % free pixels on the line became prescribed
         x_current(x_zeros)=0;   % free pixels on the line became black
         
        
     end;
 
  % STEP TWO
     
     ray_row=A(i,:);
     ray_pos=logical(ray_row~=0);
     
     x_zeros=zeros(size(free_pixels));
     x_zeros(ray_pos)=free_pixels(ray_pos);  
     x_zeros=logical(x_zeros)'; % free pixesl on the line
     
     x_attempt=x_current;
     x_attempt(x_zeros)=1; % free pixels on the line bacame white
     
     if ray_row*x_attempt==b(i),
         
         % if sum(x_attempt)>0,  disp('break');  disp(i); break; end; % control
         free_pixels(x_zeros)=0; % free pixels on the line became prescribed
         x_current(x_zeros)=1;  % free pixels on the line became white
         
         % disp('white'); disp(sum(free_hexels));
         
     end; 
     
 end;          
   i_count=i_count+1;
end;
  
x_out=x_current;

x_out(free_pixels)=0.5; % just mark free pixels

out_im=rowreshape(x_out,sqrt(n_pixels),sqrt(n_pixels))'; % out image

end
 
 
 