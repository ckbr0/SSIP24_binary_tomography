%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection_data_generator_record
% 2014
%
% Mx=c
% M -is the projection matrix
% c -is the projection vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;


%%% Load test image 
                      
            %      load ph1;   ph_orig=ph1; im_id='ph1';
            %     load ph2;   ph_orig=ph2; im_id='ph2';
            %     load ph3;   ph_orig=ph3; im_id='ph3';
            %       load ph4;   ph_orig=ph4; im_id='ph4';
            %      load ph5;   ph_orig=ph5; im_id='ph5';
            %        load ph6;   ph_orig=ph6; im_id='ph6';
            %        load ph10;   ph_orig=ph10; im_id='ph10';
            %       load ph_fat;   ph_orig=ph_fat; im_id='ph_fat';
            %       load ph_joint2;   ph_orig=ph_joint2; im_id='ph_joint2';
             %         load ph_joint1;   ph_orig=ph_joint1; im_id='ph_joint1';
            %         load ph_duck_64;   ph_orig=ph_duck_64; im_id='ph_duck_64';
            %        load ph_dog;   ph_orig=ph_dog; im_id='ph_dog';
            
            %         load ph1_4g;   ph_orig=ph1_4g; im_id='ph1_4g';
            %         load ph2_4g;   ph_orig=ph2_4g; im_id='ph2_4g';
            %         load ph3_4g;   ph_orig=ph3_4g; im_id='ph3_4g';
            %         load ph_cola;   ph_orig=ph_cola; im_id='ph_cola';
            %          load ph_hammer;   ph_orig=ph_hammer; im_id='ph_hammer';         
            %          load ph_boxer;   ph_orig=ph_boxer; im_id='ph_boxer';         
            %          load ph_cravate;   ph_orig=ph_cravate; im_id='ph_cravate'         
                      load ph_vine;   ph_orig=ph_vine; im_id='ph_vine';
                      
            %          load ph_watch;   ph_orig=ph_watch; im_id='ph_watch';
            
            
            %         load ph_elipse;   ph_orig=ph_elipse; im_id='ph_elipse';
            %         load ph_elipse_n1;   ph_orig=ph_elipse_n1; im_id='ph_elipse_n1';
            %         load ph_elipse_n2;   ph_orig=ph_elipse_n2; im_id='ph_elipse_n2';
            %          load ph_elipse_n3;   ph_orig=ph_elipse_n3; im_id='ph_elipse_n3';
            %          load ph_elipse_n4;   ph_orig=ph_elipse_n4; im_id='ph_elipse_n4';
                      
                       
            
           
               
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%% Size (n x n), must be square.
[n n]=size(ph_orig);
N=n*n; % N is the number of pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_proj=2;  % Number of projections (1,2,3,5,6 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch no_proj 
    

case 1
%%%%%%%%%%% 1 projection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_directions=135;  % in degrees, range: [0, 180]
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

proj_directions=num2str(proj_directions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
case 2
%%%%%%%%%%% 2 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_direction1=0; proj_direction2=90;


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

proj_directions=strcat(num2str(proj_direction1),'_', num2str(proj_direction2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 3
%%%%%%%%%% 3 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_directions='0_45_90';

%%%% Projection with angle 0
[A, b]= projectionmatrix(ph_orig, 0, [n/2,n/2] , n);


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


%%%%% Projection with angle pi/2
[A, b]= projectionmatrix(ph_orig, pi/2, [n/2,n/2] , n);

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

%%%%% Projection with angle pi/4
[A, b]= projectionmatrix(ph_orig, pi/4, [n/2,n/2] , n);
 
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

case 5
%%%%%%%% 5 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_directions='0_22.5_45_76.5_90';

%%%% Projection with angle 0
[A, b]= projectionmatrix(ph_orig, 0, [n/2,n/2] , n/2);
b=b(1:n); A=A(1:n,:); % just remove the quation 0=0
A=single(A);
M=A; c=b;

%%%%% Projection with angle pi/2
[A, b]= projectionmatrix(ph_orig, pi/2, [n/2,n/2] , n/2);
b=b(1:n); A=A(1:n,:); % just remove the equation 0=0
A=single(A);
M=[M;A]; c=[c;b];

%%%%% Projection with angle pi/4
[A, b]= projectionmatrix(ph_orig, pi/4, [n/2,n/2] , n/2);
A=single(A);
M=[M;A]; c=[c;b];

%%%%% Projection with angle pi/8
 [A, b]= projectionmatrix(ph_orig, pi/8, [n/2,n/2] , n/2);
 A=single(A);
 M=[M;A]; c=[c;b];

%%%%% Projection with angle 3pi/8
 [A, b]= projectionmatrix(ph_orig, 3*(pi/8), [n/2,n/2] , n/2);
 A=single(A);
 M=[M;A]; c=[c;b];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


case 6
%%%%%%%%%%%% 6 projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

proj_directions='0_30_60_90_120_150';

%%%% Projection with angle 0
[A, b]= projectionmatrix(ph_orig, 0, [n/2,n/2] , n/2);
b=b(1:n); A=A(1:n,:); % just remove the quation 0=0
A=single(A);
M =A; c=b;

%%%%% Projection with angle pi/2
[A, b]= projectionmatrix(ph_orig, pi/2, [n/2,n/2] , n/2);
b=b(1:n); A=A(1:n,:); % just remove the equation 0=0
A=single(A);
M=[M;A]; c=[c;b];

%%%%% Projection with angle pi/6
[A, b]= projectionmatrix(ph_orig, pi/6, [n/2,n/2] , n/2);
A=single(A);
M=[M;A]; c=[c;b];

%%%%% Projection with angle pi/3
[A, b]= projectionmatrix(ph_orig, pi/3, [n/2,n/2] , n/2);
A=single(A);
M=[M;A]; c=[c;b];

%%%%% Projection with angle -pi/6 (5pi/6)
[A, b]= projectionmatrix(ph_orig, -pi/6, [n/2,n/2] , n/2);
A=single(A);
M=[M;A]; c=[c;b]; 
 
%%%%% Projection with angle -pi/3 (2pi/3)
[A, b]= projectionmatrix(ph_orig, -pi/3, [n/2,n/2] , n/2);
A=single(A);
M=[M;A]; c=[c;b];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 otherwise
 disp('Invalid projection number.');
 
 end;

 
 % Make projection data record

 filename=strcat(im_id,'_p',proj_directions);

 save (filename,'M','c','im_id','ph_orig','proj_directions');
 
 
  



