function [HGN] = hologn_encoder(GNs,d,IV)
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



%An array for shifted initialization HD-vectors
E=[]; 

%For every GN shift initialization HD-vectors in the value activated in the GN
for i=1:numGN
E(end+1,:)=circshift(IV(i,:), [0 GNs(1,i)]);   
end

%Create HGN representation by majority sum on E
HGN=majority_sum(E);


end

