function [HGN, LD] = letters_encoding(d,IV)
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
%       LD low dimensional original encoding
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>


%

%Provides set of images of letters. Fig. 4 in the original paper;
load Letters 



%Initializes an array for represented images
HGN=zeros(26,d);
LD=zeros(26,35);

%Encodes every image of letter
for j=1:26 % number of letters

%Takes jth image from Letters    
pict=Letters{j,1};

%Reshapes image into pattern
pattern(1,:)=pict(:)';

%Store low dimensional representation for a pattern 
LD(j,:)=pattern;

%Creates distributed representation for a pattern across graph neurons
HGN(j,:)=hologn_encoder(pattern,d,IV);



end


end

