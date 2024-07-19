function T=im_translate(S,A,B);

% Transalte the image S by the vector A->B
% A=[x0,y0]; B=[x1,y1]

[n_row n_col]=size(S);

trans_vect=B-A;

T=zeros(n_row,n_col);

for i=1:n_row,
    for j=1:n_col,
       
        if (i+trans_vect(1)>=1) && (i+trans_vect(1)<=n_row) && (j+trans_vect(2)>=1) && (j+trans_vect(2)<=n_col),
        T(i+trans_vect(1),j+trans_vect(2))=S(i,j);
        end;
    
    end;
end;



%im_show(S);
%im_show(T);