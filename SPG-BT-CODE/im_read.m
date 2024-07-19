% Tibor Lukic, 2015

 I=imread('r_ph_cravate_sa.tif');

% I=I(50:301,50:301);

%  I=imresize(I,2.65);

%  I=I(:,10:end);
  
%  I=I(:,1:end-9);
  


% I_gray=rgb2gray(I);
 
% I_gray=I(:,:,3);

% I_gray=I_gray(1:152,1:152);


  
 I=double(I);
  
 I=I/(max(max(I)));
 
 % I=1-I;
 
 I=im_tr(I,0.5); 
  

 im_show(I);
 
 
 
 r_ph_cravate_sa=I;
 
 save r_ph_cravate_sa r_ph_cravate_sa;
 