function ConvHull=convex_hull(S)

% orientation of the set/image S in degrees
% Tibor Lukic 2017, Szeged



t=graythresh(S);
bw=im2bw(S, t);

ConvHull=bwconvhull(bw,'union');

% ConvHull=bwconvhull(bw,'objects')









