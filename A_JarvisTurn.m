function [turn]=JarvisTurn(a,b,c)

%% Function called by function Jarvis.m 
%% implementing Jarvis' march algorithm in a plane

%% INPUT
%  a, b and c are 2-dimensional column vectors.

%% OUTPUT
%  For a vector u=b-a, determine if turning to vector v=c-a is 
%  
%  clockwise    : turn = -1 = turning right
%  or 
%  anticlockwise: turn =  1 = turning left
%   
%  if the vector are colinear in opposite directions: turn = -2
%  if the vector are colinear in the same direction : turn =  2
%  if the vector are colinear and identical         : turn =  0

%% METHOD
%  Compute the sign of the sinus of the angle \theta
%  {\em measured from u to v}. Since the cross product 
%  between the two vectors u and v in a plane defines an unit vector 
%  orthogonal to u and v times the factor |u|*|v|*sin(\theta),  
%  this is equivalent to assess the sign of the determinant 
%  of the matrix formed by the vectors u and v.
    
u=b-a; v=c-a;

d=det([u v]);

turn=sign(d);

if (abs(d) < 1.e-12)        % Colinearity
    sp=u(1)*v(1)+u(2)*v(2); % Scalar product
    if (abs(sp) < 1.e-12) 
        turn = 0;           % Identical points
    else
        turn = 2*sign(sp);  % 2 * sign of scalar product  
    end 
end
    
end % End of function JarvisTurn



