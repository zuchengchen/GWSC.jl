function matrix = Rz(a)
%
% rotates coordinates of a vector around z-axis (passive rotation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix = [ cos(a),  sin(a),  0;
          -sin(a),  cos(a),  0;
             0,       0,     1 ];
