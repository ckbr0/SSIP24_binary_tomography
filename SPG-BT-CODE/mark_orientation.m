function mark_orientation(S);

% mark center of gravity and orientation
% S is an (binary) image in matrix form

[xc yc]=cent_grav(S);

alpha=orientation_vec(S(:));

xc=round(xc);
yc=round(yc);

alpha=alpha*pi/180;
k=tan(pi-alpha);
n=yc-k*xc;

[nr_row nr_col]=size(S);

pixel_density=0.4;

for x=1:pixel_density:nr_row,
    y=round(k*x+n);
    if not(y<1 || y>nr_row),
    S(y,round(x))=0.7;
    end;
end;

S(yc,xc)=0.4; S(yc-1,xc)=0.4; S(yc+1,xc)=0.4; S(yc,xc+1)=0.4; S(yc,xc-1)=0.4;
% S(yc+1,xc+1)=0.5; S(yc-1,xc+1)=0.5; S(yc+1,xc-1)=0.5; S(yc-1,xc-1)=0.5;



im_show(S);