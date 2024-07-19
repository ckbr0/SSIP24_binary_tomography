function y = row(x)

% x - matrix type: (m x n)
% y - vector type: (m*n x 1)

tmp=x.';y=tmp(:);
