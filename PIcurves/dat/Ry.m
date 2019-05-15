function matrix = Ry(a)
%
% rotates coordinates of a vector around y-axis (passive rotation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [ cos(a),  0, -sin(a);
             0,     1     0;
           sin(a),  0,  cos(a)];
