function [Res] = majority_sum(array)
% Performs bundling of array of HD-vectors
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



%Number of HD-vectors to bundle
MAX=length(array(:,1));

%Dimensionality of HD-vectors
DIM=length(array(1,:));

%If only one vector, do nothing
if MAX==1
    Res=array;
else
    
%If number of vectors is even, break ties by adding new HD-vector
if mod(MAX,2)==0
rng('default')
rng('shuffle');    
array(end+1,:)=round(rand(1,DIM));

%Increment number of HD-vectors
MAX=MAX+1;
end

%Perform majority sum
SUM=sum(array);
Res=double((SUM>MAX/2));
end

end

