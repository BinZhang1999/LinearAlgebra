function [p, e, projectionMatrix] = myProjectGF2(b, A)
% project column vector b onto the column space of A over GF(2)
% Input:
%   b: a column vector
%   A: a Matrix 
% Ouput:
%   p: the projection vector
%   e: the error vector which is perpendicular to the colume space of A
%   projectionMatrix: projectionMatrix * b = p
%
% remarks:
%   The key idea to find the projection is to find the error vector e.
%   Image a light shining from the point b to the space C(A). 
%   Light always go through the shortest way which is perpendicular to the
%   space C(A).
%   Once we find the shortest way, we find the projection.

% e = b - A*xHat
% A' * e = 0 
% A'*(b-A*xHat) = 0
% A'*b = (A'*A)*xHat
% xHat = inv(A'*A) * (A'*b)
Ax = mod(A'*A, 2);
invA = mod(myInvGF2(Ax), 2);
xHat = mod(invA * (A'*b) , 2);
e = b - A*xHat;
p = A*xHat;
% p = A*inv(A'*A)*(A'*b)
% p = (A*inv(A'*A)*A') * b
projectionMatrix = mod(A*invA*A', 2);
end