function MO=cmomt(S,p,q,cx,cy)

% centralised moment of the set/image S
% Tibor Lukic, 2015

% input: S - graysacale image in matrix form
% input: p,q - moment orders
% output: MO - moment of S

n=size(S,1); % m - number of rows
v=1:n; 

M_mat=(((v-cy).^q)' * ((v-cx).^p)).*S;

MO=sum(sum(M_mat));



