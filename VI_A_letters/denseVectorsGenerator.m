function [IV] = denseVectorsGenerator(numGN,d)
% Creates distributed representation for a pattern across graph neurons 
% 
%
% SYNOPSIS
%   HGN = sparseVectorsGenerator(numGN)
%
% DESCRIPTION
%   Creates distributed representation for a pattern across graph neurons 
%   
%   Input:
%       numGN  Number of GNs
%              
%
%   Output:
%       HGN holographic representation of the activated elements   
%           
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>

%

%Initialize an array of initialization high-dimensional vectors for every GN
%Note that one of parameters is seed, i.e. if the seed remains the same, then the
%generated HD-vectors will be the same as well
%IV=randint(numGN,d,[0,1],1);

%rng('default'); %Every time this function generates the same vectors comment the line if want to test for different vectors
IV=round(rand(numGN,d,'single')); % alterantive




end

