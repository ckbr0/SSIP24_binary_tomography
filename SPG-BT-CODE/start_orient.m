%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deterministic Binary Tomography 
% Reconstruction Method with 
% Orientation and Smothness Priors
%
% Spectral Projected Gradient / SPG type method
% 
% Tibor Lukic, 2016.
% prepared for revised version for PRL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
clear; % Clear all variables
close;

%  matlabpool; % have to be activate to reach full speed 

   
   %%%%%% START: CALCULATE PROJECTION DATA %%%%%%%%%%%%%%%%%%
   
     % load ph_dog; ph_orig=ph_dog; im_id='ph_dog';
     % load ph_joint1; ph_orig=ph_joint1; im_id='ph_joint1';
     % load ph4; ph_orig=ph4; im_id='ph4';
     % load ph6; ph_orig=ph6; im_id='ph6';
     % load ph10; ph_orig=ph10; im_id='ph10';
     % load ph10; ph_orig=ph10; im_id='ph10';
     % load ph_duck_64; ph_orig=ph_duck_64; im_id='ph_duck_64';
     %  load ph_cola; ph_orig=ph_cola; im_id='ph_cola'; % PH7
     % load ph_watch; ph_orig=ph_watch; im_id='ph_watch'; % PH8
      load ph_vine; ph_orig=ph_vine; im_id='ph_vine'; % PH9
     %  load ph_cravate; ph_orig=ph_cravate; im_id='ph_cravate'; % PH10
       
     % angles_deg=[0]; proj_directions='0'; [M c]=proj_data_calc(ph_orig,angles_deg);
       angles_deg=[0 45 73]; proj_directions='45'; [M c]=proj_data_calc(ph_orig,angles_deg);
     %  angles_deg=[90]; proj_directions='90'; [M c]=proj_data_calc(ph_orig,angles_deg);
     % angles_deg=[135]; proj_directions='135'; [M c]=proj_data_calc(ph_orig,angles_deg);
   
   %%%%%% END: CALCULATE PROJECTION DATA %%%%%%%%%%%%%%%%%%
            
disp(sprintf('Image: %s ',im_id ));
disp(sprintf('Proj. direction angles: %s ', proj_directions )); 
%%%%%%%%% END: LOAD PROJECTION DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% START: parameters of the ENERGY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%
w_proj=0.1;               % Projection weight, suggested value: 0.1, in previous versions w_proj=const=0.5
w_hom=0.5;                % Homogenity weight, suggested value: 0.5
w_orient=0;             % Orientation weight, suggested values: 0.1 or 0
%%%%%%%% END: parameters of the ENERGY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% START: parameters of the RECONSTRUCTION method %%%%%%%%%%%%%%%%%%%%%%%
tol_orient=0;             % Orientation difference toleration parameter in PERCENTAGES 
tot_it_step=5;            % Numerical gradient calculation step size
display_progress=1;       % progress visualisation, 0 or 1 
mu=0.0001;                % Initial value of the binarization parameter
mu_delta=0.05; % 0.2
Eout=10^(-3);
Ein=10^(-2); %(-2)
%Ein=10^(-7);
%%%%%%%%% END: parameters of the RECONSTRUCTION method %%%%%%%%%%%

%%%%%% START: Parameters of the SPG (inner) Algorithm %%%%%
gamma=10^(-4);   % PROBE, proposed value gamma=10^(-4)
sigma_1=0.1;
sigma_2=0.9;
m=1;             % line search parameter, m=1,2,3,4...
max_SPG_it=500;
max_tot_it=5000;
alpha_min=10^(-3);
alpha_max=10^(3);
% alpha(1)=1;
%%%%%%% END: Parameters of the SPG (inner) Algorithm %%%%%%%%%%%%

%%%% START: Initial settings %%%%%%%%%%%%%%%%%%%
[nr_proj N]=size(M); n=sqrt(N);

fprintf('Image size: %d\n',n );
fprintf('Total number of proj. rays: %d\n',nr_proj );

global orient_given;
orient_given=orientation_vec(ph_orig(:)); % global variable
time_vec=zeros(1,3);

  InitIm=0.5*ones(n,n);    % INITIAL IMAGE
% InitIm=zeros(n,n); InitIm(round(3*n/6):round(4*(n/6)),:)=0.1; InitIm=imrotate(InitIm,orient_given,'crop'); % PROBE

Xnew=InitIm(:);          % vector representation
E=ones(n*n,1);
%%%% END: Initial settings %%%%%%%%%%%%%%%%%%%%%%


if display_progress,
    im_show(ph_orig);
    set(gcf, 'Unit', 'inches'); 
    set(gcf, 'Position', [3 1.8 5 5]); % figure position and size
    set(gca, 'Unit', 'inches'); 
    set(gca, 'Position', [0.2 0.2 4.5 4.5]); % image position and size
    title(['ORIGINAL orientation: ' num2str(orientation_vec(ph_orig(:))) ]);
    figure; % settings for solution window  
    set(gcf, 'Unit', 'inches'); 
    set(gcf, 'Position', [8.5 1.8 5 5]); % figure position and size
    set(gca, 'Unit', 'inches'); 
    set(gca, 'Position', [0.2 0.2 4.5 4.5]); % image position and size
end;

wpstr=num2str(w_proj);
whstr=num2str(w_hom);
wostr=num2str(w_orient);
wtostr=num2str(tol_orient);

wsettings=strcat(' WP:',wpstr,' WH:',whstr,' WO:',wostr,' TO:',wtostr);
disp(sprintf('W_settings: %s ',wsettings ));

w_orient_given=w_orient;
       
if display_progress,   
     imshow(rowreshape(Xnew,n,n)',[0 1],'InitialMagnification','fit');
     title(['RECONSTRUCTION orientation.: ' num2str(orientation_vec(Xnew)) ' k=' num2str(0)  ]);
     drawnow; 
end;

%%%%%%%%%%% START: RECONSTRUCTION METHOD %%%%%%%%%%%%%%%%%%%%%%%%%

start_time=cputime; tic; % elapsed time measure
ng_flag=0;
total_iterations=0; 

while (total_iterations<50) || ( max(min(Xnew,1-Xnew))>Eout && norm(Xold-Xnew,1)>0 && total_iterations<max_tot_it), % GENERAL STOPPING CRITERIUM
    
    % clear fun_bt_values; 
    fun_bt_values=zeros(max_SPG_it+1);
    
    orient_x=orientation_vec(Xnew);
    if isnan(orient_x), 
          w_orient=0;
        else
          w_orient=w_orient_given;
    end;
    
   
    fun_bt_values(1)=energy(Xnew,n,n,M,c,w_proj,w_hom,w_orient,mu);
        
    k=1; alpha_k=1; 
    
    %%%%%%%% START: gradient calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g=gradient_hom_term_vec(Xnew,n,n);   
    if not(isnan(orient_x)) && (w_orient~=0),
        if ( (abs(orient_x-orient_given)/orient_given)*100>tol_orient && mod(total_iterations,tot_it_step)==0 ) || (ng_flag==0),
        ng_orient=num_gradient_par('orientation_vec',Xnew); ng_flag=1;
        end;   
        gg=w_proj*M'*(M*Xnew-c) + w_hom*g + w_orient*(orient_x-orient_given)*ng_orient - mu*(Xnew-0.5*E); 
    else 
        gg=w_proj*M'*(M*Xnew-c) + w_hom*g  - mu*(Xnew-0.5*E); 
    end;
    %%%%%%%% END: gradient calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%% START: SPG ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%
    
    while  (k<2) || (  norm( proj(Xnew-gg)-Xnew,inf) > Ein &&  k < max_SPG_it && lambda>10^(-6) ), 
        
    
    
    %%%%%%%%%%%%%%% new direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prj=proj(Xnew-alpha_k*gg);
    dk=proj(Xnew-alpha_k*gg)-Xnew; % Scaled Projected Gradient
    %%%%%%%%%%%%%%%  new deriction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%% START: LINE SEARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    fun_max=max( fun_bt_values( 1+k-min(k,m) : k ) );
    
    x_plus=Xnew+dk;    
    delta=dk'*gg;
    lambda=1;  % PROBE lambda=1
    
    while energy(x_plus,n,n,M,c,w_proj,w_hom,w_orient,mu) > fun_max+gamma*lambda*delta,
        
              
        lambda_temp=-0.5*lambda^2*delta/(energy(x_plus,n,n,M,c,w_proj,w_hom,w_orient,mu)-fun_bt_values(k)-lambda*delta );
        if and( (lambda_temp >= sigma_1) , (lambda_temp<= sigma_2*lambda) ),
            lambda=lambda_temp;
        else
            lambda=0.5*lambda;
        end;
            x_plus=Xnew+lambda*dk; % NEW ATEMPT
            
        % fprintf('LEFT:  %f',  energy_bt_orient(x_plus,n,n,M,c,w_hom,w_orient,mu)  ); % PROBE
        % fprintf('  RIGHT:  %f',  fun_max+gamma*lambda*delta  ); % PROBE 
        % fprintf('  <dk,gg>:  %f', delta   ); % PROBE
        % fprintf('  norm(gg):  %f', norm(gg)   ); % PROBE
        % fprintf('  norm(dk):  %f', norm(dk)   ); % PROBE
        % fprintf('  fun_max:  %f\n',  fun_max  ); % PROBE 
        % fprintf('  lambda:  %f\n',  lambda  ); % PROBE
             
    end;
    
    %%%%%%%%%%%%%%%%%% END: LINE SEARCH  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    Xold=Xnew;
    Xnew=x_plus; % NEW ITERATION
    
    k=k+1; total_iterations=total_iterations+1; 
    
    if display_progress,   
     imshow(rowreshape(Xnew,n,n)',[],'InitialMagnification','fit');
     title(['RECONSTRUCTION orientation.: ' num2str(orientation_vec(Xnew)) ' k=' num2str(k-1) ' it=' num2str(total_iterations) ]);
     drawnow; 
    % for zz=1:1000000, fprintf('delay counter:  %d \n' ,zz ); end; %PROBE 
    end;
    
    orient_x=orientation_vec(Xnew);
    if isnan(orient_x), 
        w_orient=0;
    else
        w_orient=w_orient_given;
    end;
    
    fun_bt_values(k)=energy(Xnew,n,n,M,c,w_proj,w_hom,w_orient,mu);
    s_k=Xnew-Xold;
     
    %%%%%%%% START: gradient calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g=gradient_hom_term_vec(Xnew,n,n); 
    if not(isnan(orient_x)) && (w_orient~=0),
        if ( (abs(orient_x-orient_given)/orient_given)*100>tol_orient && mod(total_iterations,tot_it_step)==0) || (ng_flag==0),
        ng_orient=num_gradient_par('orientation_vec',Xnew); ng_flag=1;        
        end;
        gg_new=w_proj*M'*(M*Xnew-c) + w_hom*g + w_orient*(orient_x-orient_given)*ng_orient - mu*(Xnew-0.5*E); 
    else
        gg_new=w_proj*M'*(M*Xnew-c) + w_hom*g - mu*(Xnew-0.5*E);
    end;
    %%%%%%%% END: gradient calculation %%%%%%%%%%%%%%%%%%%%%%%
      
    y_k=gg_new - gg; gg=gg_new;
    beta_k=s_k'*y_k;
    
    if (beta_k <= 0),
        alpha_k=alpha_max;
    else
        alpha_k=min(alpha_max,max(alpha_min,(s_k'*s_k)/beta_k));
    end;
     
    end;
    
    %%%%%%%%%%% END: SPG ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%

%   w_hom=w_hom*0.4; % PROBE
    
    mu=mu+mu_delta; % increasing the binarization parameter
%   mu=2*mu+mu_delta;
%   mu=1.05*mu;

   
   
   if display_progress,
      elapsed_time=cputime-start_time;
      time_vec(1)=fix(elapsed_time/3600);
      time_vec(2)=fix((elapsed_time-time_vec(1)*3600)/60);
      time_vec(3)=fix((elapsed_time-time_vec(2)*60-time_vec(1)*3600));
      fprintf('binarization term:  %f\n' , mu); 
      fprintf('total and inner it.:  %d %d\n' , total_iterations,k-1); 
      fprintf('change:  %f\n' , norm(Xold-Xnew,1) ); 
      fprintf('elapsed time: %dh %dm %ds\n  %f' , time_vec); 
      disp(' ');
   end;

end;

   out_im=rowreshape(Xnew,n,n)'; % Shape it back to an image
   out_im=im_tr(out_im,0.5);     % Tresholding
   
   
   

%%%%%%%%%%% END: RECONSTRUCTION METHOD %%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% START: DISPLAY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

elapsed_time=toc;
time_vec(1)=fix((elapsed_time)/3600);
time_vec(2)=fix((elapsed_time-time_vec(1)*3600)/60);
time_vec(3)=fix((elapsed_time-time_vec(2)*60-time_vec(1)*3600));
projection_error=norm(M*out_im(:)-c);
pixel_error=sum(sum(abs(ph_orig-out_im)));
rnmp=pixel_error/sum(sum(abs(ph_orig)));
MR=(pixel_error/N)*100;
orient_orig=orientation_vec(ph_orig(:));
orient_rec=orientation_vec(out_im(:));
orient_diff=abs(orient_rec-orient_orig);
cent_grav_diff=norm(cent_grav(ph_orig)-cent_grav(out_im));



imshow(out_im,[],'InitialMagnification','fit');
title(['RECONSTRUCTION orientation.:' num2str(orient_rec) '  it=' num2str(total_iterations) ]);
drawnow;

[cpe cmr cent_orient_diff trans_im]=cpe_calc(ph_orig, out_im);

disp(sprintf('Image name:  %s ',im_id ));
disp(sprintf('Projection direction angles:  %s ', proj_directions ));
fprintf('Orientation difference toleration:    %d\n', tol_orient );
fprintf('Line search parameter M:    %d\n', m ); 
fprintf('w_hom: %f\n', w_hom);
fprintf('w_orient: %f\n', w_orient);
fprintf('Total number of iterations: %d\n', total_iterations );
fprintf('Elapsed time: %dh %dm %ds\n' , time_vec);
fprintf('Projection error PRE:           %f\n',  projection_error  );
fprintf('Pixel error PE:       %d pixels\n',  pixel_error   );
fprintf('Relative number of misclassified pixels, rNMP:       %f \n', rnmp    );
fprintf('Misclassified rate (in percentage), MR:       %f \n', MR    );
fprintf('\n');
fprintf('ERRORS MR/PRE:  %f/%f \n',  MR,projection_error  );
fprintf('Orient_diff: %f degrees\n', orient_diff );
fprintf('\n');
fprintf('Orientation of the original: %f degrees\n', orient_orig );
fprintf('Orientation of the reconstruction: %f degrees\n', orient_rec );
fprintf('Orient_diff: %f degrees\n', orient_diff );
fprintf('Center of gravity difference: %f \n', cent_grav_diff );
fprintf('\n');
fprintf('Centered pixel errorors CPE/CMR:  %d/%f \n',  cpe,cmr   );
fprintf('Centered orientation difference:    %f \n',  cent_orient_diff  );

%if display_progress,
%    imshow(out_im,[],'InitialMagnification','fit');
%    title('RECONSTRUCTION');
%    drawnow;     
%else
%    im_show(ph_orig);
%    set(gcf, 'Unit', 'inches'); 
%    set(gcf, 'Position', [2 1.8 5 5]); % figure position and size
%    set(gca, 'Unit', 'inches'); 
%    set(gca, 'Position', [0.2 0.2 4.5 4.5]); % image position and size
%    title('ORIGINAL');
     
%    im_show(out_im);
%    set(gcf, 'Unit', 'inches'); 
%    set(gcf, 'Position', [8 1.8 5 5]); % figure position and size
%    set(gca, 'Unit', 'inches'); 
%    set(gca, 'Position', [0.2 0.2 4.5 4.5]); % image position and size
%    title('RECONSTRUCTION');    
%end;

%%%%%%%%%%%%%%%%%% END: DISPLAY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%% START: SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_im=abs(ph_orig-out_im);
for i=1:n,
    for j=1:n,
        if diff_im(i,j)==1,
            if ph_orig(i,j)==1,
                diff_im(i,j)=0.5;
            end;
        end;
    end;
end;
    

if w_orient~=0,
    filename=strcat('r_',im_id,'_p',proj_directions,'_WP',wpstr,'_WH',whstr,'_WO',wostr,'_TO',wtostr,'.mat');
else
    filename=strcat('r_',im_id,'_p',proj_directions,'_WP',wpstr,'_WH',whstr,'_WO',wostr,'.mat');
end;

save (filename,'out_im','ph_orig','nr_proj','diff_im','proj_directions','projection_error','pixel_error','rnmp',...
               'MR','cpe','cmr','time_vec','total_iterations', 'orient_rec', 'orient_orig', 'orient_diff', 'display_progress','tol_orient',...
               'tot_it_step', 'w_proj','w_hom', 'w_orient' );
           
%%%%%%%%%%%%%%%%%%% END: SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           


% mark_orientation(ph_orig);
% mark_orientation(out_im);

    
