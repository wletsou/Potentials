% find the psuedo-inverse of the n-choose-2 by n matrix B of combinations of
% pairs

function out = combination_matrix(n,varargin)

p = inputParser;
p.CaseSensitive = 1;

parse(p);

Pairs = nchoosek(1:n,2);
A = eye(n);

B = (1 / (n - 2) ) * cell2mat(transpose(arrayfun(@(I)...
    sum(A(setdiff(1:n,Pairs(I,:)),:),1),1:nchoosek(n,2),'UniformOutput',false)));...
    % position of ones in indicates the numbers included in the combination

BTBinv = (eye(n) * (-( (n - 2) * (n - 3) + (n - 1) ) - (n - 3)) + ones(n) * (n - 3) ) /...
    ( -2 * (n - 2) ) * 2 * (n - 2) / (n - 1); % inverse of B^T * B

out = BTBinv * transpose(B);