% CENTERED PIXEL ERROR
% inputs: out_im, ph_orig

% PH1 / ph_dog 
 
 
  %load ph_dog_p0;
  % load r_ph_dog_p0_O-;
 %  load r_ph_dog_p0_O+;

  %load ph_dog_p45;
  %load r_ph_dog_p45_O-;
  %load r_ph_dog_p45_O+;

  %load ph_dog_p90;
  %load r_ph_dog_p90_O-;
 % load r_ph_dog_p90_O+;

 % load ph_dog_p135;
  %load r_ph_dog_p135_O-;
 % load r_ph_dog_p135_O+;

% PH2 / ph_joint1

 % load ph_joint1_p0;
  %load r_ph_joint1_p0_O-;
 % load r_ph_joint1_p0_O+;

  % load ph_joint1_p45;
  % load r_ph_joint1_p45_O-;
  %load r_ph_joint1_p45_O+;

 % load ph_joint1_p90;
  %load r_ph_joint1_p90_O-;
 % load r_ph_joint1_p90_O+;

 % load ph_joint1_p135;
  %load r_ph_joint1_p135_O-;
 % load r_ph_joint1_p135_O+;


% PH3/ ph4

   load ph4_p0;
  % load r_ph4_p0_O-;
   load r_ph4_p0_O+;

% load ph4_p45;
% load r_ph4_p45_O-;
%load r_ph4_p45_O+;

 % load ph4_p90;
  %load r_ph4_p90_O-;
 % load r_ph4_p90_O+;

% load ph4_p135;
 %load r_ph4_p135_O-;
% load r_ph4_p135_O+;

% PH4/ ph6

% load ph6_p0;
 %load r_ph6_p0_O-;
% load r_ph6_p0_O+;

 %load ph6_p45;
 %load r_ph6_p45_O-;
 %load r_ph6_p45_O+;

 %load ph6_p90;
 %load r_ph6_p90_O-;
 %load r_ph6_p90_O+;

 %load ph6_p135;
 %load r_ph6_p135_O-;
 %load r_ph6_p135_O+;

% PH5/ ph10

  %load ph10_p0;
  %load r_ph10_p0_O-;
  %load r_ph10_p0_O+;

%  load ph10_p45;
 %load r_ph10_p45_O-;
% load r_ph10_p45_O+;

 % load ph10_p90;
 % load r_ph10_p90_O-;
 %load r_ph10_p90_O+;

 %load ph10_p135;
 %load r_ph10_p135_O-;
 %load r_ph10_p135_O+;

% PH6/ ph_duck_64

 %load ph_duck_64_p0;
 %load r_ph_duck_64_p0_O-;
 %load r_ph_duck_64_p0_O+;

 %load ph_duck_64_p45;
 %load r_ph_duck_64_p45_O-;
 %load r_ph_duck_64_p45_O+;

%load ph_duck_64_p90;
%load r_ph_duck_64_p90_O-;
%load r_ph_duck_64_p90_O+;

% load ph_duck_64_p135;
 %load r_ph_duck_64_p135_O-;
% load r_ph_duck_64_p135_O+;




[x_o y_o]=cent_grav(ph_orig);
[x_r y_r]=cent_grav(out_im);
x_o=round(x_o); y_o=round(y_o);
x_r=round(x_r); y_r=round(y_r);
T=im_translate(out_im, [y_r x_r], [y_o x_o]);


projection_error=norm(M*T(:)-c);
pixel_error=sum(sum(abs(ph_orig-T)));
rnmp=pixel_error/sum(sum(abs(ph_orig)));
orient_orig=orientation(ph_orig);
orient_rec=orientation(T);
orient_diff=abs(orient_rec-orient_orig);

diff_im=abs(ph_orig-T);
n=size(T);
for i=1:n,
    for j=1:n,
        if diff_im(i,j)==1,
            if ph_orig(i,j)==1,
                diff_im(i,j)=0.5;
            end;
        end;
    end;
end;

% DISPLAY RESULTS
fprintf('\n');
disp(sprintf('Image name:  %s ',im_id ));
disp(sprintf('Projection direction angles:  %s ', proj_directions ));



fprintf('\n');

fprintf('PE/rNM/PRE      %d  / %f / %f \n',  pixel_error, rnmp*100, projection_error);
fprintf('Orient_diff: %f degrees\n', orient_diff );

fprintf('\n');

%  mark_orientation(ph_orig); mark_orientation(out_im); mark_orientation(T); % im_show(diff_im);



% END: DISPLAY RESULTS



