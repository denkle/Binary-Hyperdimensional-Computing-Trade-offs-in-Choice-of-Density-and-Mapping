function [HGN] = sparse_letters_encoding(IV,d)
% Creates distributed representation (HGN) for every image of letter  
% 
%
% SYNOPSIS
%   HGN=letters_encoding(  )
%
% DESCRIPTION
%   Creates distributed representation (HGN) for every image of letter  
%   
%   Input:
%       There are no input parameters      
%
%   Output:
%       HGN array of distributed representations for every image of letter
%           


%

%Provides set of images of letters. Fig. 4 in the original paper;
load Letters 

%Set the dimensionality of HD-vectors
%d=2048;

%Initializes an array for represented images
HGN=single(zeros(26,d));

%Encodes every image of letter
for j=1:26 % number of letters

%Takes jth image from Letters    
pict=Letters{j,1};

%Reshapes image into pattern
pattern(1,:)=pict(:)';

%Creates distributed representation for a pattern across graph neurons
HGN(j,:)=sparse_hologn_encoder_IV(pattern,IV);

end


end

