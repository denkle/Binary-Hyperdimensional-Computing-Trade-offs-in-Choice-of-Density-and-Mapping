function [ decoded, nz ] = item_memory_c(vector, memory )
% Recall the closest vector in item-memory and outputs corresponding index
%
% SYNOPSIS
%   decoded = item_memory_c(vector, memory)
%
% DESCRIPTION
%   Recall the closest vector in item-memory and outputs corresponding
%   value.
%
%   Input:
%       vector  vector to be recalled
%       memory  item memory of HoloGN       
%       Please note that inputs are assumed to be arrays of complex numbers with
%       possible values 0 or j
%
%   Output:
%       decoded index in memory, which HD-vector has the smallest Hamming   
%           distance with the input vector
% 
% AUTHOR
%   Denis Kleyko <denis.kleyko@ltu.se>


    % Hamming distances between the input and the values in item-memory
    HD=imag(memory*vector);
    
    % Index of the closest vector in item memory
    [v]=min(HD);
    nz=nnz(HD==v); %number of nonzero elements
    
    if nz==1
        [~,decoded]=min(HD);
    else
        
        V = find(HD==v);
        rng('shuffle');
        ind=randi([1,nz],1,1);
        decoded=V(ind);        
    end

   
end

