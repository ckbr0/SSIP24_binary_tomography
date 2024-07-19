function [M, c]=proj_data_calc(ph_orig,angles_deg)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection_data_generator_record
% 2016, Tibor Lukic.
% 
%
% Mx=c
% M -is the projection matrix
% c -is the projection vector
% angles_deg - angles of the projection directions in degrees [0, 180] (vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Size (n x n), must be square.
[n n]=size(ph_orig);
N=n*n; % N is the number of pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_proj=size(angles_deg,2); % number of projections 



switch no_proj 
    

case 1
%%%%%%%%%%% 1 projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_directions=angles_deg;  % in degrees, range: [0, 180]
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);  % probe n/2

rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);

A=single(A);
M =A; c=b;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
case 2
%%%%%%%%%%% 2 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_direction1=angles_deg(1); proj_direction2=angles_deg(2);
% in degrees, range: [0, 180]

%%%%%%%%%% proj 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proj_directions=proj_direction1;  % in degrees, range: [0, 180]
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);  % probe n/2

rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);
M =A; c=b;
%%%%%%%%%% end: proj 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% proj 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proj_directions=proj_direction2;  % in degrees, range: [0, 180]
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);  % probe n/2

rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);
%%%%%%%%%% end: proj 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=[M;A]; c=[c;b];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 3
%%%%%%%%%% 3 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in degrees, range: [0, 180]
proj_direction1=angles_deg(1); 
proj_direction2=angles_deg(2);
proj_direction3=angles_deg(3);


%%%% Projection 1
proj_directions=proj_direction1;  
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);

rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);
M =A; c=b;


%%%%% Projection 2
proj_directions=proj_direction2;  
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);


rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);

M=[M;A]; c=[c;b];

%%%%% Projection 3
proj_directions=proj_direction3;  
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);

 
rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);
M=[M;A]; c=[c;b];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


otherwise
    
M= zeros(no_proj*n,n*n,'single');


    
proj_direction1=angles_deg(1); 

%%%% Projection 1
proj_directions=proj_direction1;  
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);

rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);
nr_rowA=size(A,1);
M(1:nr_rowA,:)=A; c=b;
nr_rowM=nr_rowA;
%%%% END Projection 1
    
for pr_count=2:no_proj,
    
%%%%% Projections
proj_directions=angles_deg(pr_count);
if proj_directions>90, 
    direction_angle=proj_directions-180; 
else
   direction_angle=proj_directions;
end;
direction_angle=(direction_angle*pi)/180; % convert to radians

clear A; clear b;
[A, b]= projectionmatrix(ph_orig, direction_angle, [n/2,n/2] , n);


rowA=size(A,1);  % just remove equations 0=0
nonzero_projections=ones(rowA,1);    
for i=1:rowA,
    if A(i,:)==0, nonzero_projections(i)=0;
    end;      
end;
nonzero_projections=logical(nonzero_projections);
A=A(nonzero_projections,:);
b=b(nonzero_projections);
A=single(A);

nr_rowA=size(A,1);

M(nr_rowM+1:nr_rowM+nr_rowA,:)=A;

nr_rowM=nr_rowM+nr_rowA;

c=[c;b];
            
%%% END Projections


end;
     
       
 
 end;

 
 
 
 
  



