function [ array_c ] = bin2com(array)
% Transforms the array or vector of binary values into array or vector or
% complex values
%
% SYNOPSIS
%   array_c  = bin2com(array)
%
% DESCRIPTION
%   Transforms the array or vector of binary values into array or vector or
%   complex values
%
%   Input:
%       array  array of binary numbers     
%
%   Output:
%       array_c array of complex numbers   
%           
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>
%

    % Create an inversion of array
    array_inv=abs(array-1);
    
    % Create an array of complex numbers where original 1 corresponds to 1 
    % and original 0 corresponds to j. New array has single precision
    array_c=single(complex(array,array_inv));
    

end

