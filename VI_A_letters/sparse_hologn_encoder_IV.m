function [HGN] = sparse_hologn_encoder_IV(GNs,IV)
% Creates distributed representation for a pattern across graph neurons 
% 
%
% SYNOPSIS
%   HGN = hologn_encoder(GNs)
%
% DESCRIPTION
%   Creates distributed representation for a pattern across graph neurons 
%   
%   Input:
%       GNs  pattern to be represented by HoloGN
%              
%
%   Output:
%       HGN holographic representation of the activated elements   
%           
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>


%

%Number of GNs is determined by the legnth of the input pattern
numGN=size(GNs,2);

%Initialize an array of initialization sparse high-dimensional vectors for every GN
%IV = sparseVectorsGenerator(numGN);


%An array for shifted initialization HD-vectors
E=[]; 

%For every GN shift initialization HD-vectors in the value activated in the GN
for i=1:numGN
E(end+1,:)=circshift(IV(i,:), [0 GNs(1,i)]);   
end

%Create sHoloGN representation by CDT on E
%CHECK IF NEED CDTWhile
%M=300; % Approximate number of ones in the distributed representation of pattern
%HGN=CDTWhile(E,M);
HGN=single(CDT(E));
end

