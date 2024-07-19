function [cx,cy]=cent_grav(S)

% center of gravity of the set/image S
% Tibor Lukic 2012

% input: S - graysacale image in matrix form
% output: cx,cy - coordinates of the center of gravity

mom00=momt(S,0,0);

cx=momt(S,1,0)/mom00;
cy=momt(S,0,1)/mom00;