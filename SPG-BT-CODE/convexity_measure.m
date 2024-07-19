function cm=convexity_measure(S)

% orientation of the set/image S in degrees
% Tibor Lukic 2017, Szeged



t=graythresh(S);
bw=im2bw(S, t);

ConvHull=bwconvhull(bw,'union');

ConvHullArea = bwarea(ConvHull);
bwHullArea = bwarea(bw);

cm=bwHullArea/ConvHullArea;











