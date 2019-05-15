function matrix = Rx(a)
%
% rotates coordinates of a vector around x-axis (passive rotation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [ 1,    0,      0;
           0,  cos(a), sin(a);
           0, -sin(a), cos(a) ];
