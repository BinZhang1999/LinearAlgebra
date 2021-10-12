function invA = myInvGF2(A)
% find the inv matrix over GF(2)
%
% Input: 
%   A: the full rank square matrix over GF(2)
% Ouput:
%   invA: inv matrix of A

[mRow, nCol] = size(A);
invA = nan(mRow);
% if a square and full rank matrix
if mRow ~= nCol
    disp('Not a square matrix !');
    return; 
end
if rank(gf(A)) ~= mRow
    disp('Not a full rank matrix !');
end
% perform invert
augMatrix = nan(mRow, 2*mRow);
augMatrix(:,1:mRow) = A; augMatrix(:,1+mRow:end) = eye(mRow); 
[matrixRowEchelon, indexColPivot, rankOfMatrix] = ...
    myEliminateGF2(augMatrix);
matrixEliminated =  myBackSubstitutionGF2(matrixRowEchelon, ...
    indexColPivot, rankOfMatrix);
invA = matrixEliminated(:,1+mRow:end);

end