function display_err(ph_orig,out_im)


N=size(ph_orig,1)*size(ph_orig,2);


pixel_error=sum(sum(abs(ph_orig-out_im)));
rnmp=pixel_error/sum(sum(abs(ph_orig)));
MR=(pixel_error/N)*100;
orient_orig=orientation_vec(ph_orig(:));
orient_rec=orientation_vec(out_im(:));
orient_diff=abs(orient_rec-orient_orig);
cent_grav_diff=norm(cent_grav(ph_orig)-cent_grav(out_im));

[cpe cmr cent_orient_diff trans_im]=cpe_calc(ph_orig, out_im);


fprintf('Pixel error PE:       %d pixels\n',  pixel_error   );
fprintf('Relative number of misclassified pixels, rNMP:       %f \n', rnmp    );
fprintf('Misclassified rate (in percentage), MR:       %f \n', MR    );
fprintf('\n');

fprintf('Orient_diff: %f degrees\n', orient_diff );
fprintf('\n');
fprintf('Orientation of the original: %f degrees\n', orient_orig );
fprintf('Orientation of the reconstruction: %f degrees\n', orient_rec );
fprintf('Orient_diff: %f degrees\n', orient_diff );
fprintf('Center of gravity difference: %f \n', cent_grav_diff );
fprintf('\n');
fprintf('Centered pixel errorors CPE/CMR:  %d/%f \n',  cpe,cmr   );
fprintf('Centered orientation difference:    %f \n',  cent_orient_diff  );
