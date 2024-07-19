function ori=orientation_vec(v)

% orientation of the set/image in degrees
% Tibor Lukic 2015

% input: v - graysacale image in vector form
% input: n - sqrt of the length of v <=> number of rows and columns 
% output: ori - orientation of v, degrees in range [0 180]

n=sqrt(length(v));
S=rowreshape(v,n,n)'; % Shape it back to an image

[cx cy]=cent_grav(S);

cm11=cmomt(S,1,1,cx,cy);
cm20=cmomt(S,2,0,cx,cy);
cm02=cmomt(S,0,2,cx,cy);

% mom=2*cm11/(cm20-cm02);

or1=(1/2)*atan(2*cm11/(cm20-cm02));
or2=or1+pi/2;

F1=sin(or1)^2*cm20+cos(or1)^2*cm02-sin(2*or1)*cm11;
F2=sin(or2)^2*cm20+cos(or2)^2*cm02-sin(2*or2)*cm11;

if (F1<=F2), 
    ori=or1; 
else
    ori=or2; 
end;

ori=ori*(180/pi);

if ori>0, ori=180-ori;
else
    ori=-ori;
end;


