function diff_im=difference_im(ph_orig,out_im)

n=size(ph_orig,1);
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