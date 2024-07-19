function [cpe cmr orient_diff trans_im]=cpe_calc(ph_orig, out_im);

% CENTERED PIXEL ERROR
% inputs: out_im, ph_orig

% output: 
% cpe= centered pixel error
% cmr= centered misclassification rate (in percentage)
% trans_im= translated image

% Tibor Lukic, 2016.


[x_o y_o]=cent_grav(ph_orig);
[x_r y_r]=cent_grav(out_im);
x_o=round(x_o); y_o=round(y_o);
x_r=round(x_r); y_r=round(y_r);
trans_im=im_translate(out_im, [y_r x_r], [y_o x_o]);



cpe=sum(sum(abs(ph_orig-trans_im)));
N=size(ph_orig,1)*size(ph_orig,1);
cmr=(cpe/N)*100;

orient_orig=orientation(ph_orig);
orient_rec=orientation(trans_im);
orient_diff=abs(orient_rec-orient_orig);







