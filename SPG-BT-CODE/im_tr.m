function S=im_tr(I,v)

S=I;

S(S>v)=1;
S(S<=v)=0;




