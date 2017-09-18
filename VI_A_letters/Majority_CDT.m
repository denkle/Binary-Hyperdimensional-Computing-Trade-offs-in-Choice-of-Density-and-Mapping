function [Res] = Majority_CDT(array,M_in,d)
% Performs Context-Dependent Thinning operation
% 
%
% SYNOPSIS
%   Res = maj_n(array)
%
% DESCRIPTION
%   Performs bundling of array of HD-vectors
%   
%
%   Input:
%       array  array of HD-vector for bundling; rows represent HD-vectors;
%              values in array are 0s ans 1s.
%
%   Output:
%       Res resultant HD-vector   
%           
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>
%



%Number of HD-vectors to bundle
MAX=length(array(:,1));

%Dimensionality of HD-vectors
DIM=length(array(1,:));

%Efficient implementation of disjunctive superposotion under the whole array
superposition=sum(array,1);

Mexp=d*(1-(1-M_in/d)^35)^2;



%set the density
if sum(double(superposition>=1))<Mexp
    M=sum(double(superposition>=1));
else
    M=Mexp; %number of ones in the vector
end 


Res=zeros(1,DIM); % resultant vector
for i=1:M
    [v,decoded]=max(superposition);
    %disp([v,decoded])
    Res(1, decoded)=1;
    superposition(1, decoded)=-1;
end


end

