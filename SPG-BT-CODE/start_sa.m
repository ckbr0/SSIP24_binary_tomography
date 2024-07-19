% SIMULATED ANNEALING algorithm
% for Binary Tomography
% 
% A Benchmark Evaluation of Large-Scale OPtimization Approaches to Binary Tomography
% S. Weber, A. Nagy and T. Schule




%%%%%% START: CALCULATE PROJECTION DATA %%%%%%%%%%%%%%%%%%
   
     % load ph_dog; ph_orig=ph_dog; im_id='ph_dog';
     % load ph_joint1; ph_orig=ph_joint1; im_id='ph_joint1';
     % load ph4; ph_orig=ph4; im_id='ph4';
     % load ph6; ph_orig=ph6; im_id='ph6';
     % load ph10; ph_orig=ph10; im_id='ph10';
     % load ph10; ph_orig=ph10; im_id='ph10';
     % load ph_duck_64; ph_orig=ph_duck_64; im_id='ph_duck_64';
     %  load ph_cola; ph_orig=ph_cola; im_id='ph_cola'; % PH7
        load ph_watch; ph_orig=ph_watch; im_id='ph_watch'; % PH8
     % load ph_vine; ph_orig=ph_vine; im_id='ph_vine'; % PH9
     %  load ph_cravate; ph_orig=ph_cravate; im_id='ph_cravate'; % PH10
       
       angles_deg=[0]; proj_directions='0'; [M c]=proj_data_calc(ph_orig,angles_deg);
     %  angles_deg=[45]; proj_directions='45'; [M c]=proj_data_calc(ph_orig,angles_deg);
     %  angles_deg=[90]; proj_directions='90'; [M c]=proj_data_calc(ph_orig,angles_deg);
     % angles_deg=[135]; proj_directions='135'; [M c]=proj_data_calc(ph_orig,angles_deg);
     %  angles_deg=[0 90]; proj_directions='0 90'; [M c]=proj_data_calc(ph_orig,angles_deg);
   
   %%%%%% END: CALCULATE PROJECTION DATA %%%%%%%%%%%%%%%%%%

    
   start_im=brute_force(M,c);
   
   start_im_vec=start_im(:);
   free_pixels=logical(start_im_vec==0.5);
   nr_free_pixels=sum(free_pixels);
   
   fprintf('SIMULATED ANNEALING ALGORITHM  \n' );
   
     
%%%%%%%% START: parameters of the ENERGY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%
w_proj=0.1;               % Projection weight, suggested value: 0.1, in previous versions w_proj=const=0.5
w_hom=5;                % Homogenity weight, suggested value: 0.5
%%%%%%%% END: parameters of the ENERGY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%




 % PARAMETERS OF THE SA ALGORITHM
 T_start=4;
 T_min=10^(-10);
 T_factor=0.9;
 R_objective=0.00001;
 
 nr_restart=1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 N=size(M,2); n=sqrt(N);
 
 start_im_vec(start_im_vec==0.5)=0;
 InitIm=start_im_vec; 
 
 x_old=InitIm; x_new=x_old;
 
 C_start=energy_sa(x_old,n,n,M,c,w_proj,w_hom);  
 C_old=C_start;
 
 
 figure; % settings for solution window  
 set(gcf, 'Unit', 'inches'); 
 set(gcf, 'Position', [8.5 1.8 5 5]); % figure position and size
 set(gca, 'Unit', 'inches'); 
 set(gca, 'Position', [0.2 0.2 4.5 4.5]); % image position and size
 
 T=T_start;
 
 for restart_count=1:nr_restart,
     
 
 C_start=energy_sa(x_old,n,n,M,c,w_proj,w_hom);  
 C_old=C_start;
 
 while (T >= T_min) && (C_old/C_start >= R_objective)
     
    for i=1:2*N,
       
       x_new=x_old;  
        
       ia=randi(nr_free_pixels,1); % discrete uniform distribution on 1:nr_free_pixels    
       xcfh=x_new(free_pixels); 
       xcfh(ia)=1-xcfh(ia);        % change the color of the random pixel
       x_new(free_pixels)=xcfh;    % new image attempt
     
        
       C_new=energy_sa(x_new,n,n,M,c,w_proj,w_hom); 
       z=rand; % random number from [0,1], uniform distribution
       delta_C=C_new-C_old;
       if delta_C<0 || exp(-delta_C/T)>z,
           x_old=x_new; % accept changes
           C_old=C_new;
       end;
       
       
     end; % end for
     
     T=T*T_factor;
     fprintf('T= %2.16f \n', T );
     
     imshow(rowreshape(x_old,n,n)',[0 1],'InitialMagnification','fit');
     title([ im_id '  RESTART=' num2str(restart_count) ]);
     drawnow;
     
 end;
 
 T=0.5;
 
 % T=T_start/(restart_count+1);
 
 end;
 
 
 out_im=rowreshape(x_old,n,n)';
 
 display_errors(ph_orig,out_im,M,c);
 
 im_show(out_im);
 
 diff_im=difference_im(ph_orig,out_im);
 
 im_show(diff_im);

 
 
 filename=strcat('r_',im_id,'_p',proj_directions,'.mat');
 
 save (filename,'out_im');
 
 
 
 
 