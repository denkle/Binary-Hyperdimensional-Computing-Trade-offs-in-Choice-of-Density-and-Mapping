function [IV] = sparseVectorsGenerator_sparsity(numGN,d,rho)
% Creates initialization vectors for the given dimensionality and sparsity 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>

%
%d - sets the dimensionality of HD-vectors
%rho - expected density of HD-vector

%Initialize an array of initialization high-dimensional vectors for every GN

%rng('default'); %Every time this function generates the same vectors comment the line if want to test for different vectors
rng('shuffle')
IV=rand(numGN,d);
IV=single(IV<=rho);
end

