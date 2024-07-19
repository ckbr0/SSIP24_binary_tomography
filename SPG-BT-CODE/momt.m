function MO=momt(S,p,q)

% moment of the set/image S
% Tibor Lukic, 2012

% input: S - graysacale image in matrix form
% input: p,q - moment orders
% output: MO - moment of S

n=size(S,1); % m - number of rows, n - number of collumns
v=1:n;

M_mat=((v.^q)' * (v.^p)).*S;

MO=sum(sum(M_mat));



