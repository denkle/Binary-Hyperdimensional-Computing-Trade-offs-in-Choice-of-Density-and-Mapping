function [ decoded ] = item_memory_ed(vector, memory )
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


    % Euclidean distances between the input and the values in item-memory
    %HD=(memory*vector);
    vector=vector';   
    HD=zeros(26,1);
    for i=1:26
       HD(i,1)=sqrt(sum((memory(i,:)-vector).^2));
    end
    
    
    % Index of the closest vector in item memory
    [v]=min(HD);
    nz=nnz(HD==v); %number of nonzero elements
    
    if nz==1
        [~,decoded]=min(HD);
    else
        V = find(HD==v);
        rng('shuffle')
        decoded=V(randi([1,nz],1,1));        
    end
    
    

end

