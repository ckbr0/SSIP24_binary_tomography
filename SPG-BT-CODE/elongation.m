function EL=elongation(S)

% elongation of the set/image S
% Tibor Lukic 2012

% input: S - graysacale image in matrix form
% output: EL - elongation of S


max_elong=0.5*(cmomt(S,2,0)+cmomt(S,0,2)+sqrt(4*cmomt(S,1,1)^2+(cmomt(S,2,0)-cmomt(S,0,2))^2));
min_elong=0.5*(cmomt(S,2,0)+cmomt(S,0,2)-sqrt(4*cmomt(S,1,1)^2+(cmomt(S,2,0)-cmomt(S,0,2))^2));

% constrains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if min_elong<10^(-3), min_elong=10^(-3); end; 
  if max_elong>10^(3), max_elong=10^(3); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EL=max_elong/min_elong;

% EL=10*log(abs(EL)); 