% Novi Sad, 2016.

extended_p=312; % new total number of projection rays;

n_proj=size(M,1);

proj_count=1;

A=M;

while size(A,1)<extended_p,
 
 b=M(proj_count,:);
 A=[b;A];
 
 proj_count=proj_count+2;
    
end;

M=A;

c=M*ph_orig(:);

clear extended_p proj_count n_proj A b;