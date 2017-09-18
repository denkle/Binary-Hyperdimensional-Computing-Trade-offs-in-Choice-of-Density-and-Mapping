function message = holoGN
  assignin('base','lookupItemMemory', @lookupItemMemory);
  assignin('base','genRandomHV', @genRandomHV);
  assignin('base','cosAngle', @cosAngle);
  assignin('base','computeNgram_4N_features', @computeNgram_4N_features);
  assignin('base','computeNgram_1N_feature_perm', @computeNgram_1N_feature_perm);
  assignin('base','hdctrain', @hdctrain); 
  assignin('base','initItemMemories', @initItemMemories);
  assignin('base','genTrainData', @genTrainData);
  assignin('base','downSampling', @downSampling);
  assignin('base','hdcpredict', @hdcpredict);
  
  assignin('base','test_slicing', @test_slicing);
  assignin('base','hdcpredict_window', @hdcpredict_window);
  assignin('base','majority_sum', @majority_sum);
  
    %functions related to sparse encoding
  assignin('base','initItemMemoriesSparse', @initItemMemoriesSparse);
  assignin('base','computeNgram_4N_featuresSparse', @computeNgram_4N_featuresSparse);
  assignin('base','computeCDT', @computeCDT);
  assignin('base','hdctrainSparse', @hdctrainSparse);
  assignin('base','MajorityCDT', @MajorityCDT);
  assignin('base','hdcpredictSparse', @hdcpredictSparse);
     %END functions related to sparse encoding
     %function to initialize iM for HoloGN and other encodings (does not matter sparse or dense)
  assignin('base','initItemMemoriesHoloGN', @initItemMemoriesHoloGN);
  assignin('base','initItemMemoriesNonLinear', @initItemMemoriesNonLinear);
  assignin('base','initItemMemoriesLinear', @initItemMemoriesLinear);
  assignin('base','initItemMemoriesProjection', @initItemMemoriesProjection);
  assignin('base','hdctrainProjection', @hdctrainProjection);
  assignin('base','hdcpredictProjection', @hdcpredictProjection);
  
   %Support functions for running conventinal machine learningn methods
  assignin('base','featurizeData', @featurizeData); 
  assignin('base','featurizeTestingData', @featurizeTestingData); 
  assignin('base','make_centorid', @make_centorid); 
  
  message='Done importing functions to workspace'; 
end




function [L_SAMPL_DATA, SAMPL_DATA] = genTrainData (data, labels, trainingFrac, order)
%
% DESCRIPTION   : generate a dataset useful for training using a fraction of input data 
%
% INPUTS:
%   data        : input data
%   labels      : input labels
%   trainingFrac: the fraction of data we shouls use to output a training dataset
%   order       : whether preserve the order of inputs (inorder) or randomly select
%   donwSampRate: the rate or stride of downsampling
% OUTPUTS:
%   SAMPL_DATA  : dataset for training
%   L_SAMPL_DATA: corresponding labels
%    

	rng('default');
    rng(1);
    L1 = find (labels == 1);
    L2 = find (labels == 2);
    L3 = find (labels == 3);
    L4 = find (labels == 4);
    L5 = find (labels == 5);
	L6 = find (labels == 6);
	L7 = find (labels == 7);
   
    L1 = L1 (1 : floor(length(L1) * trainingFrac));
    L2 = L2 (1 : floor(length(L2) * trainingFrac));
    L3 = L3 (1 : floor(length(L3) * trainingFrac));
    L4 = L4 (1 : floor(length(L4) * trainingFrac));
    L5 = L5 (1 : floor(length(L5) * trainingFrac));
	L6 = L6 (1 : floor(length(L6) * trainingFrac));
	L7 = L7 (1 : floor(length(L7) * trainingFrac));
 
    if order == 'inorder'
		Inx1 = 1:1:length(L1);
		Inx2 = 1:1:length(L2);
		Inx3 = 1:1:length(L3);
		Inx4 = 1:1:length(L4);
		Inx5 = 1:1:length(L5);
		Inx6 = 1:1:length(L6);
		Inx7 = 1:1:length(L7);
	else
		Inx1 = randperm (length(L1));
		Inx2 = randperm (length(L2));
		Inx3 = randperm (length(L3));
		Inx4 = randperm (length(L4));
		Inx5 = randperm (length(L5));
		Inx6 = randperm (length(L6));
		Inx7 = randperm (length(L7));
    end
    
    
    L_SAMPL_DATA = labels (L1(Inx1));
    L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L2(Inx2)))];
    L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L3(Inx3)))];
    L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L4(Inx4)))];
    L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L5(Inx5)))];
	L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L6(Inx6)))];
	L_SAMPL_DATA = [L_SAMPL_DATA; (labels(L7(Inx7)))];
    %L_SAMPL_DATA = L_SAMPL_DATA';
    
    SAMPL_DATA   = data (L1(Inx1), :);
    SAMPL_DATA   = [SAMPL_DATA; (data(L2(Inx2), :))];
    SAMPL_DATA   = [SAMPL_DATA; (data(L3(Inx3), :))];
    SAMPL_DATA   = [SAMPL_DATA; (data(L4(Inx4), :))];
    SAMPL_DATA   = [SAMPL_DATA; (data(L5(Inx5), :))];
	SAMPL_DATA   = [SAMPL_DATA; (data(L6(Inx6), :))];
	SAMPL_DATA   = [SAMPL_DATA; (data(L7(Inx7), :))];
end

function [downSampledData, downSampledLabels] = downSampling (data, labels, donwSampRate)
%
% DESCRIPTION   : apply a downsampling to get rid of redundancy in signals 
%
% INPUTS:
%   data        : input data
%   labels      : input labels
%   donwSampRate: the rate or stride of downsampling
% OUTPUTS:
%   downSampledData: downsampled data
%   downSampledLabels: downsampled labels
%    
	j = 1;
	for i = 1:donwSampRate:length(data)
		downSampledData (j,:) = data(i, :);
		downSampledLabels (j) = labels(i);
        j = j + 1;
    end
    
    downSampledLabels = downSampledLabels';
end
	
function randomHV = genRandomHV(D,density)
%
% DESCRIPTION   : generate a random vector with zero mean 
%
% INPUTS:
%   D           : Dimension of vectors
% OUTPUTS:
%   randomHV    : generated random vector

    if mod(D,2)
        disp ('Dimension is odd!!');
    else
        %randomIndex = randperm (D);
        %randomHV (randomIndex(1 : D/2)) = 1;
        %global codingMode;
	    %switch codingMode
        %    case 'dense_bipolar'
        %        randomHV (randomIndex(D/2+1 : D)) = -1;
        %    case 'dense_binary'
        %        randomHV (randomIndex(D/2+1 : D)) = 0;
        %end
        
        %randomHV=round(rand(1,D));
        randomHV=double(rand(1,D)>=(1-density));
        
        global codingMode;
	    switch codingMode
            case 'dense_bipolar'     
                randomHV=2*(randomHV-0.5);
        end
        
        
        
    end
end

function [iM] = initItemMemories (D, MAXL, N,seed,density)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
% seed          : seed to set random generator in a prespecified position. This will randomize results due to initial vectors
% density       : set initial density of 1s for vectors.
% OUTPUTS:
%   iM          : item memory
 
	iM = containers.Map ('KeyType','char','ValueType','any');
    rng ('default');
    rng (seed);
    
    % Here we form N*4 features, hence we require N*4 orthogonal initial vectors        
    for j = 1:1:N*4
	    currentHV = genRandomHV (D,density);
	    %randomIndex = randperm (D);
        %DK. New part to implement for varying densities
        density_act=sum(currentHV==1); %actual number of ones
        pos1=find(currentHV==1);
        pos0=find(currentHV~=1);
        randomIndex1 = randperm (density_act);
        randomIndex0 = randperm (D-density_act);
        
        % Iterate over dicrete levels and 'continuously' generate vectors for related values	
        for i = 0:1:MAXL
            %disp(i)
            key = strcat(int2str(i),'x',int2str(j));
            iM(key) = currentHV;
            
		    %D / 2 / MAXL = 238
            %SP = floor(D/2/(MAXL)); %Added -1. Before the dot product betwee number 0 and number 21 was a bit more than 0
            %DK. New part to implement for varying densities
            SP=  floor((density_act-D*((density_act/D)^2))/MAXL);
            
            % start index and end index for flipping bits
		    startInx = (i*SP) + 1;
		    endInx = ((i+1)*SP); %DK. There was a small typo. I have deleted +1
            % flip these random bits
            
            if i~=MAXL
            global codingMode;
	        switch codingMode
                case 'dense_bipolar'
		            %currentHV (randomIndex(startInx : endInx)) = currentHV (randomIndex(startInx : endInx)) * -1;
                    currentHV (pos1(randomIndex1(startInx : endInx))) = -1;
                    currentHV (pos0(randomIndex0(startInx : endInx))) = 1;
                    
                    
                case 'dense_binary'
		            %currentHV (randomIndex(startInx : endInx)) = not (currentHV (randomIndex(startInx : endInx)));
                     currentHV (pos1(randomIndex1(startInx : endInx))) = 0;
                    currentHV (pos0(randomIndex0(startInx : endInx))) = 1;
                    
                    
            end
            end
        end
    end
end


function randomHV = lookupItemMemory (itemMemory, channelValue, channelID, percision)
%
% DESCRIPTION   : recalls a vector from item Memory (here continous) based on inputs
%
% INPUTS:
%   itemMemory  : item memory
%   channelValue: the analog value of EMG channel
%   channelID   : ID of EMG channel from 1 to 4*N
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   randomHV    : return the related vector

    % Discretization of EMG signal based on required percision
    quantized = int64 (channelValue * percision);
    % generates a key to find related vector for quantized value in channel with this ID
    key = strcat (int2str(quantized), 'x', int2str(channelID));
    if itemMemory.isKey (key) 
        randomHV = (itemMemory (key));
        %fprintf ('READING KEY: %s\n', key);
    else
        fprintf ('CANNOT FIND THIS KEY: %s\n', key);       
    end
end

function ngram = computeNgram_4N_features (block, iM, D, N, percision)
%
% DESCRIPTION   : compute an ngram vector for a block of input data; a block has 4*N sampels data
%               This encoding uses 4*N separate features and hence 4*N item memories
%
% INPUTS:
%   block       : input data
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   ngram       : the ngram vector
    
    ngram =( zeros (4*N,D));
            for t = 1:1:N
            % Go over the 4 channels of EMG
                for c = 1:1:4
			        %ngram = ngram + lookupItemMemory (iM, block(t,c), 4*(t-1)+c, percision);
                     ngram(4*(t-1)+c,:)=lookupItemMemory (iM, block(t,c), 4*(t-1)+c, percision); 
                end
            end
    
    global bipolarDistretizationMode
    global codingMode;
	switch codingMode
        case 'dense_bipolar'           
             if strcmp(bipolarDistretizationMode,'on')
                ngram=(majority_sum(ngram));
             else
                ngram=sum(ngram);
             end
        case 'dense_binary'
            %ngram=sum(ngram);
            ngram=(majority_sum(ngram));            
    end
end

 
% function ngram = computeNgram_1N_feature_perm (block, iM, D, N, percision)
% %
% % DESCRIPTION   : compute an ngram vector for a block of input data; a block has 4*N sampels data
% %               This encoding uses 4 separate features and hence 4 item memories and then applies permutation
% % INPUTS:
% %   block       : input data
% %   iM          : item memory
% %   D           : Dimension of vectors
% %   N           : size of n-gram, i.e., window size
% %   percision   : percision used in quantization of input EMG signals
% %
% % OUTPUTS:
% %   ngram       : the ngram vector
%     
%     global codingMode;
% 	switch codingMode
%         case 'dense_bipolar'
%         	ngram = zeros (1,D);
%             % Go over various time stamps to capture signal history (n-gram)
%             for t = 1:1:N
%                 % Compute a sum vector for each time stamp of signal, i.e., a vertical slicing of EMG data (N=1)
% 	            sumv = zeros (1,D);
%                 for c = 1:1:4
% 			        sumv = sumv + lookupItemMemory (iM, block(t,c), c, percision);
%                 end
% 	            % Now, permure it to compute n-gram = sumv(1) + p(sumv(2)) + pp(sumv(3) + ... 
%                 ngram = ngram + circshift (sumv, [1, t-1]);
%             end
% 
%         case 'dense_binary'
%             ngram = zeros (1,D);
%             % Go over various time stamps to capture signal history (n-gram)
%             for t = 1:1:N
%                 % Compute a sum vector for each time stamp of signal, i.e., a vertical slicing of EMG data (N=1)
% 	            sumv = zeros (1,D);
%                 for c = 1:1:4
%                     retrieveVec = lookupItemMemory (iM, block(t,c), c, percision);
% 			        sumv = [sumv; retrieveVec];
%                 end
%                 sumv (1,:) = [];
%                 sumv = mode (sumv);
% 	            % Now, permure it to compute n-gram = sumv(1) + p(sumv(2)) + pp(sumv(3) + ... 
%                 rotatedVec = circshift (sumv, [1, t-1]);
%                 ngram = [ngram; rotatedVec];
%             end
%             ngram(1,:) = [];
%             ngram = mode(ngram);
%     end
% end

             
function [numPat, AM] = hdctrain (labels, data, iM, D, N, percision, cuttingAngle) 
%
% DESCRIPTION   : train an associative memory based on input training data
%
% INPUTS:
%   labels      : training labels
%   data        : EMG training data
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%   cuttingAngle: threshold angle for not including a vector into SUM vector
%
% OUTPUTS:
%   AM          : Trained associative memory
%   numPat      : Number of stored patterns for each class of AM
%

    %global vectors_training
    %vectors_training=[];
    global trainingMode;
    global bipolarDistretizationMode;

     		% initialize an empty AM
			AM = containers.Map ('KeyType','double','ValueType','any');
			%fprintf ('Total traning samples size = %d\n', length(labels));
			
			for label = 1:1:max(labels)
                
				AM (label) = zeros (1,D);
				numPat (label) = 0;
			end

			% walk through the input data 
			i = 1;
			while i < length(labels)-N+1
				% if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
				if labels(i) == labels(i+N-1)		
					% compute n-gram; the encoderMode determines the type of encoding algorithm
					global encoderMode;
					switch encoderMode
						case '4N_features'
                            
							ngram = computeNgram_4N_features (data(i : i+N-1,:), iM, D, N, percision);
						case 'N_feature_perm'
							ngram = computeNgram_1N_feature_perm (data(i : i+N-1,:), iM, D, N, percision);
					end

					
                    
						AM (labels(i+N-1)) = AM (labels(i+N-1)) + ngram;
						numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
                        
                        
                    % walk through data with stride 1
					i = i + 1;
				else
					% jump to the next window that has a 'stable' gesture where all labels are equal
					i = i + N - 1; %DK. double check
				end
            end
            
    switch trainingMode    
       case 'integer_vect'
            if strcmp(bipolarDistretizationMode,'on')
                for label = 1:1:max(labels)
                    AM (label) =majority_sum_sup(AM (label),numPat (label)) ;	 
                end     
            end
       case 'binary_vect'    
           	for label = 1:1:max(labels)
				AM (label) =majority_sum_sup(AM (label),numPat (label)) ;	 
			end
    end      
            
    
        
% 	switch trainingMode
%         case 'integer_vect'
%      		% initialize an empty AM
% 			AM = containers.Map ('KeyType','double','ValueType','any');
% 			%fprintf ('Total traning samples size = %d\n', length(labels));
% 			%DK. commented for simulation
% 			for label = 1:1:max(labels)
% 				AM (label) = zeros (1,D);
% 				numPat (label) = 0;
% 			end
% 
% 			% walk through the input data 
% 			i = 1;
% 			while i < length(labels)-N+1
% 				% if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
% 				if labels(i) == labels(i+N-1)		
% 					% compute n-gram; the encoderMode determines the type of encoding algorithm
% 					global encoderMode;
% 					switch encoderMode
% 						case '4N_features'
% 							ngram = computeNgram_4N_features (data(i : i+N-1,:), iM, D, N, percision);
% 						case 'N_feature_perm'
% 							ngram = computeNgram_1N_feature_perm (data(i : i+N-1,:), iM, D, N, percision);
% 					end
% 
% 					% check whether the produced n-gram is already in AM
% 					%DK. I have trurned off the cosine similarity check for
% 					%the sake of fair comparison. Alternatively we have to
% 					%do this check for binary vectors
%                     %angle = cosAngle(ngram, AM (labels(i+N-1)));
% 					%if angle < cuttingAngle | isnan(angle)
% 						AM (labels(i+N-1)) = AM (labels(i+N-1)) + ngram;
% 						numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
% 					%end
% 					
%                     % walk through data with stride 1
% 					i = i + 1;
% 				else
% 					% jump to the next window that has a 'stable' gesture where all labels are equal
% 					i = i + N - 1; %DK. double check
% 				end
%             end
%             
%             if strcmp(bipolarDistretizationMode,'on')
%                 for label = 1:1:max(labels)
%                     sigHV=AM(label);
%                     buf=zeros(1,D);
%                     buf=buf+(sigHV>0);
%                     buf=buf-(sigHV<0);
%                     bufi=find(buf==0);
%                     b_per=randperm(length(bufi));
%                     buf(bufi(b_per(1:round(length(bufi)/2))))=1;  
%                     buf(bufi(b_per(round(length(bufi)/2)+1:end)))=-1;                     
%                     %AM (label) = 2*(double(AM(label)>0 )  -0.5)  ;
%                     AM (label) = buf  ;
%                     
%                 end
%                 
%             end
%             
%             
% 		
% 		case 'binary_vect'
%      		% initialize an empty AM
% 			AM = containers.Map ('KeyType','double','ValueType','any');
% 			%fprintf ('Total traning samples size = %d\n', length(labels));
% 			%DK. Commented for the simulations 
% 		
% 			
% 			for label = 1:1:max(labels)
% 				%AM (label) = zeros (1,D);
% 				numPat (label) = 0;
% 				trainVecList = (zeros (1,D));
% 				% walk through the input data 
% 				i = 1;
% 				while i < length(labels)-N+1
% 					% if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
% 					if (labels(i) == labels(i+N-1)) && (labels(i) == label)		
% 						% compute n-gram; the encoderMode determines the type of encoding algorithm
% 						global encoderMode;
% 						switch encoderMode
% 							case '4N_features'
% 								ngram = computeNgram_4N_features (data(i : i+N-1,:), iM, D, N, percision);
% 							case 'N_feature_perm'
% 								ngram = computeNgram_1N_feature_perm (data(i : i+N-1,:), iM, D, N, percision);
% 						end
% 
% 						% check whether the produced n-gram is already in AM
% 						%angle = cosAngle(ngram, mode (trainVecList));
% 						%if angle < cuttingAngle | isnan(angle)
% 							trainVecList = [trainVecList; ngram]; %DK for the fair comparison we should include fultering as well or turn it off for the bipolar case
% 							numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
% 						%end
% 						% walk through data with stride 1
% 						i = i + 1;
% 					else
% 						% jump to the next window that has a 'stable' gesture where all labels are equal
% 						i = i + N - 1;
% 					end
% 				end
% 				trainVecList(1,:) = [];
% 				%AM (label) = mode (trainVecList);
%                 AM (label) =(majority_sum(trainVecList));
%                 
%                 
% 			end
%     end	
%     
    
	%DK. commented for simulations		
    %for label = 1:1:max(labels)
	%	fprintf ('Class= %d \t mean= %0.3f \t n_pat=%d \t created \n', label, mean(AM(label)), numPat(label));
    %end
end
              

function [accExcTrnz, accuracy] = hdcpredict (labelTestSet, testSet, AM, iM, D, N, percision)
%
% DESCRIPTION   : test accuracy based on input testing data
%
% INPUTS:
%   labelTestSet: testing labels
%   testSet     : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy for all situations
%   accExcTrnz  : classification accuracy excluding the transitions between gestutes
%
	correct = 0;        % count the number of correct classifications
    numTests = 0;       % count the total number of tests
	tranzError = 0;     % count the number of tests that went wrong but were between gesture transitions 
	%STAT=zeros(1,5);    %DK. To test simple hypotesis
    
    %global testngrams;
    global bipolarDistretizationMode;
    k=1;
    % Go over the test data with stride of N
	for i = 1:N:length(testSet)-N+1
		numTests = numTests + 1;
		% set actual label as the majority of labels in the window!
        actualLabel = mode(labelTestSet (i : i+N-1));
        
        % compute ngram vector for the testing window; the encoderMode determines the type of encoding algorithm   
        global encoderMode;
		switch encoderMode
            case '4N_features'
		        sigHV = computeNgram_4N_features (testSet (i : i+N-1,:), iM, D, N, percision);
            case 'N_feature_perm'
		        sigHV = computeNgram_1N_feature_perm (testSet (i : i+N-1,:), iM, D, N, percision);
        end
            
  
       global codingMode;
        switch codingMode
            case 'dense_bipolar'
                if strcmp(bipolarDistretizationMode,'on')
                maxAngle = -1;
                predicLabel = -1;
                for label = 1:1:max(labelTestSet)
                    %angle = cosAngle(AM (label), sigHV);
                    angle = sum(AM (label).* sigHV);
                    if (angle > maxAngle)
                        maxAngle = angle;
                        predicLabel = label;
                    end
                end    
                    
                else
                maxAngle = -1;
                predicLabel = -1;
                for label = 1:1:max(labelTestSet)
                    angle = cosAngle(AM (label), sigHV);
                    if (angle > maxAngle)
                        maxAngle = angle;
                        predicLabel = label;
                    end
                end
                end
                
                
            case 'dense_binary'
                maxAngle = 1;
                predicLabel = -1;
                for label = 1:1:max(labelTestSet)
                    angle = sum(xor(AM (label), sigHV))/D; %Hamming distance
                    if (angle < maxAngle)
                        maxAngle = angle;
                        predicLabel = label;
                    end
                end

        
        
        end
        
                if predicLabel == actualLabel
                    correct = correct + 1;
                    
                elseif labelTestSet (i) ~= labelTestSet(i+N-1)
                    tranzError = tranzError + 1;
                end
        
        %STAT(1,actualLabel)=STAT(1,actualLabel)+1;
    end
    %disp(STAT)
    accuracy = correct / numTests;
	accExcTrnz = (correct + tranzError) / numTests;
end


function [accuracy] = test_slicing (labels, data, AM, iM, D, N, percision)
%
% DESCRIPTION   : test accuracy based on input testing data but it automatically slices the test signal to different gestures
%
% INPUTS:
%   labels      : testing labels
%   data        : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy when using sliced windows
%

	correct = 0;
    numTests = 0;
    start = length(labels)-1;	
    
    % Go over the test data with stride of 1 and finds a window during which all labels are equal (hence, set start and stop indices) 
    for i = 1:1:length(labels)-1
         if (labels (i) == labels (i+1)) & (start > i)
            start = i;
         elseif (labels (i) ~= labels (i+1)) & (start <= i)
            stop = i;
            window = stop - start;
            window = max (window, N);
		
            if start < 1 | stop+window > length(labels)-1
        		fprintf (' !!!! selected gesture is out of range %d:%d !!! \n', start, stop);
            else
                % Find the label for this window (that might have multiple ngrams) by picking up the ngram that has highest similarity
				[maxAngle, predicLabel] = hdcpredict_window_max (data(start : start+window, :), AM, iM, D, N, percision);
                	
		        numTests = numTests + 1;
                if predicLabel == labels(start)
			        correct = correct + 1;
                else
        		    fprintf (':-( %d went to %d : for signal range of %d:%d \n', labels(start), predicLabel, start, stop);
                end
            end
            start = length(data)-1;
        end
    end

    accuracy = correct / numTests;
end

function [maxAngle, predicLabel] = hdcpredict_window_max (testSet, AM, iM, D, N, percision)
%
% DESCRIPTION   : predicts a label for an ngram that has the highest similarity within a fixed window size 
%
% INPUTS:
%   testSet     : The block of data in which this function searches to find the highest similarity
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   maxAngle    : the maximum angle of ngram found in the window
%   predicLabel : predicted label for the ngram
%

    maxAngle = -1;
    predicLabel = -1;

    % Walk through the block of data (that might have multiple ngrams) with stride 1
    for i = 1:1:length(testSet(:,1))-N+1
        % compute ngram vectors that are available in  the testing window; the encoderMode determines the type of encoding algorithm   
        global encoderMode;
		switch encoderMode
            case '4N_features'
		        sigHV = computeNgram_4N_features (testSet(i : i+N-1,:), iM, D, N, percision);
            case 'N_feature_perm'
        		sigHV = computeNgram_1N_feature_perm (testSet(i : i+N-1,:), iM, D, N, percision);
        end
        % Pick up the label for ngram that has the maximum simularity with stored patterns in AM
	    for label = 1:1:5
			angle = cosAngle (AM(label), sigHV);
			if (angle > maxAngle)
				maxAngle = angle;
				predicLabel = label;
			end
        end
    end
end


function cosAngle = cosAngle (u, v)
    cosAngle = dot(u,v)/(norm(u)*norm(v));
end


function [Res] = majority_sum(array)
% Performs bundling of array of HD-vectors
% SYNOPSIS
%   Res = maj_n(array)
% DESCRIPTION
%   Performs bundling of array of HD-vectors
%   Input:
%       array  array of HD-vector for bundling; rows represent HD-vectors;
%              values in array are 0s ans 1s.
%   Output:
%       Res resultant HD-vector  
%Number of HD-vectors to bundle


MAX=length(array(:,1));
%Dimensionality of HD-vectors
DIM=length(array(1,:));

global codingMode;
switch codingMode
case 'dense_binary'  
    
%If only one vector, do nothing
if MAX==1
    Res=array;
else 
%If number of vectors is even, break ties by adding new HD-vector
if mod(MAX,2)==0
%rng('shuffle');    
array(end+1,:)=round(rand(1,DIM));
%Increment number of HD-vectors
MAX=MAX+1;
end

%Perform majority sum
SUM=sum(array);
Res=double((SUM>MAX/2));
end



case 'dense_bipolar'     
        %TBD
        if MAX==1
         Res=array;
        else 
        %If number of vectors is even, break ties by adding new HD-vector
        if mod(MAX,2)==0
        %rng('shuffle');    
        array(end+1,:)=2*( round(rand(1,DIM))-0.5 );
        %Increment number of HD-vectors
        MAX=MAX+1;
        end

        %Perform majority sum
        SUM=sum(array);
        Res=2*(double((SUM>0))-0.5);
        end     
end
end

function [Res] = majority_sum_sup(array,MAX)
% Performs bundling of array of HD-vectors
% SYNOPSIS
%   Res = maj_n(array)
% DESCRIPTION
%   Performs bundling of array of HD-vectors
%   Input:
%       array  array of HD-vector for bundling; rows represent HD-vectors;
%              values in array are 0s ans 1s.
%   Output:
%       Res resultant HD-vector  
%Number of HD-vectors to bundle


%MAX=length(array(:,1));
%Dimensionality of HD-vectors
DIM=length(array(1,:));

global codingMode;
switch codingMode
case 'dense_binary'  
    
        %If only one vector, do nothing
        if MAX==1
            Res=array;
        else 
        %If number of vectors is even, break ties by adding new HD-vector
        if mod(MAX,2)==0
        %rng('shuffle');    
        array=array+round(rand(1,DIM));
        %Increment number of HD-vectors
        MAX=MAX+1;
        end

        %Perform majority sum
        Res=double((array>MAX/2));
        end



case 'dense_bipolar'     
        %TBD
        if MAX==1
         Res=array;
        else 
        %If number of vectors is even, break ties by adding new HD-vector
        if mod(MAX,2)==0
        %rng('shuffle');    
        array=array+ 2*( round(rand(1,DIM))-0.5 );
        %Increment number of HD-vectors
        MAX=MAX+1;
        end
        %Perform majority sum
        Res=2*(double((array>0))-0.5);
        end     
end
end




%% Sparse encoder


function [iM] = initItemMemoriesSparse (D, MAXL, N, rho,seed)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
%   rho           :density of HD-vector 

% OUTPUTS:
%   iM          : item memory
 
	iM = containers.Map ('KeyType','char','ValueType','any');
    rng ('default');
    rng (seed);
    
    % Here we form N*4 features, hence we require N*4 orthogonal initial vectors        
    for j = 1:1:N*4
	    %currentHV = genRandomHV (D);
        
        M=round(rho*D); %number of ones in a sparse HD-vector
        
        currentHV=double(rand(1,D)>=(1-rho));
        %disp([j,sum(currentHV)])
        
        density_act=sum(currentHV==1); %actual number of ones
        pos1=find(currentHV==1); %find indexes of non zero elements
        pos0=find(currentHV~=1);%find indexes of  zero elements
        randomIndex1 = randperm (density_act);
        randomIndex0 = randperm (D-density_act);
        
        
        % Iterate over dicrete levels and 'continuously' generate vectors for related values	
        for i = 0:1:MAXL
            key = strcat(int2str(i),'x',int2str(j));
            iM(key) = currentHV;
		    %D / 2 / MAXL = 238
            %SP = floor(Mr/MAXL);
            SP=  floor((density_act-D*((density_act/D)^2))/MAXL);
            % start index and end index for flipping bits
		    startInx = (i*SP) + 1;
		    endInx = ((i+1)*SP);
            % flip these random bits
		    %currentHV (randomIndex(startInx : endInx)) = currentHV (randomIndex(startInx: endInx)) * -1;
            if i~=MAXL
            %currentHV (randomIndex_nz(startInx : endInx)) = 0; %null several nonzero elements
            %currentHV (randomIndex_z(startInx : endInx)) = 1;  %activate several zero elements
            currentHV (pos1(randomIndex1(startInx : endInx))) = 0;%null several nonzero elements
            currentHV (pos0(randomIndex0(startInx : endInx))) = 1;%activate several zero elements
            end
            
        end
     
    end
end

function ngram = computeNgram_4N_featuresSparse (block, iM, D, N, percision)
%
% DESCRIPTION   : compute an ngram vector for a block of input data; a block has 4*N sampels data
%               This encoding uses 4*N separate features and hence 4*N item memories
%
% INPUTS:
%   block       : input data
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   ngram       : the ngram vector, which is already binarized and thinned (thinning is done twice)

	ngram = zeros (1,D);
    % Go over various time stamps to capture signal history (n-gram)
    for t = 1:1:N
        % Go over the 4 channels of EMG
        for c = 1:1:4
			ngram = ngram + lookupItemMemory (iM, block(t,c), 4*(t-1)+c, percision);
        end
    end
    %create binary vector
    ngram=(double(ngram>=1));
    ngram = computeCDT(ngram); 

    
end


function ngram = computeCDT(ngram)
% Performs Context-Dependent Thinning operation
% 
%   Input:
%       binary vector before thinning
%
%   Output:
%        resultant HD-vector   

%Dimensionality of HD-vectors
DIM=length(ngram(1,:));

coef_shift=round(0.8*DIM);

%Efficient implementation of disjunctive superposotion under the whole array
superposition=ngram;

%Implementation of fixed number of permutations for thinning
ngram=zeros(1,DIM); % resultant vector
K=2;% number of permutations, i.e. thinning operations
for i=1:K
    shifted=circshift(superposition, [0 coef_shift+i]);
    thinned=and(superposition,shifted);
    ngram=or(ngram,thinned); %update of the resultant vector    
    %disp([sum(thinned),sum(Res)]);
end

ngram=double(ngram);

end

function [Res] = MajorityCDT(superposition,M)
% SYNOPSIS
%   Res = maj_n(array)
% DESCRIPTION
%   Performs binarization of superposition of sparse HD-vectors
%   Input:      superposition  superposition of sparse HD-vectors
%   Output:      Res resultant binary sparse HD-vector   

%Dimensionality of HD-vectors
DIM=length(superposition(1,:));

%set the density
if sum(double(superposition~=0))<M
    M=sum(double(superposition~=0));
end 

Res=zeros(1,DIM); % resultant vector
[~,decoded]=sort(superposition,'descend');
Res(1, decoded(1:M))=1;
end



function [numPat, AM_b] = hdctrainSparse (labels, data, iM, D, N, percision, cuttingAngle) 
%
% DESCRIPTION   : train an associative memory based on input training data
%
% INPUTS:
%   labels      : training labels
%   data        : EMG training data
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%   cuttingAngle: threshold angle for not including a vector into SUM vector
%
% OUTPUTS:
%   AM          : Trained associative memory
%   numPat      : Number of stored patterns for each class of AM
%
    
    d_ng=0;% counter density of ngrams. will be used in majorutyCDT
    c_ng=0;%counter for the number of ngrams. will be used in majoruty sum
    
    % initialize an empty AM
    AM = containers.Map ('KeyType','double','ValueType','any');
    AM_b = containers.Map ('KeyType','double','ValueType','any'); %binarized AM
	%fprintf ('Total traning samples size = %d\n', length(labels));
    for label = 1:1:max(labels)
    	AM (label) = zeros (1,D);
        AM_b(label) = zeros (1,D);
	    numPat (label) = 0;
    end

    % walk through the input data 
    i = 1;
    while i < length(labels)-N+1
        % if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
        if labels(i) == labels(i+N-1)		
            % compute n-gram; the encoderMode determines the type of encoding algorithm
           % global encoderMode;
	        %switch encoderMode
                %case '4N_features'
			        ngram = computeNgram_4N_featuresSparse (data(i : i+N-1,:), iM, D, N, percision);
               % case 'N_feature_perm'
			        %ngram = computeNgram_1N_feature_perm (data(i : i+N-1,:), iM, D, N, percision);
            %end

			% check whether the produced n-gram is already in AM
            %angle = cosAngle(ngram, AM (labels(i+N-1)));
            %angle = sum(AM_b (labels(i+N-1)).* ngram)/sum(ngram);

			%if angle < cuttingAngle | isnan(angle)
				AM (labels(i+N-1)) = AM (labels(i+N-1)) + ngram;
                numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
                d_ng=d_ng+sum(ngram);
                c_ng=c_ng+1;
                d_mean=round(d_ng/c_ng);
                %AM_b (labels(i+N-1)) =MajorityCDT( AM (labels(i+N-1)),d_mean);
                	            
			%end
            % walk through data with stride 1
            i = i + 1;
		else
            % jump to the next window that has a 'stable' gesture where all labels are equal
            i = i + N - 1;
            %i=i+1;
        end
    end
        d_mean=round(d_ng/c_ng);
        
    for label = 1:1:max(labels)
             AM_b (label) =MajorityCDT( AM (label),round(0.9*d_mean));
            %disp(round(0.9*d_mean));
		%fprintf ('Class= %d \t mean= %.0f \t n_pat=%d \t created \n', label, mean(AM(label)), numPat(label));
    end
end

function [accExcTrnz, accuracy] = hdcpredictSparse (labelTestSet, testSet, AM, iM, D, N, percision)
%
% DESCRIPTION   : test accuracy based on input testing data
%
% INPUTS:
%   labelTestSet: testing labels
%   testSet     : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy for all situations
%   accExcTrnz  : classification accuracy excluding the transitions between gestutes
%
	correct = 0;        % count the number of correct classifications
    numTests = 0;       % count the total number of tests
	tranzError = 0;     % count the number of tests that went wrong but were between gesture transitions 
	
    %numpat=zeros(1,5);
    %numpatpred=zeros(1,5);
    
    % Go over the test data with stride of N
	for i = 1:N:length(testSet)-N+1
		numTests = numTests + 1;
		% set actual label as the majority of labels in the window!
        actualLabel = mode(labelTestSet (i : i+N-1));
        
        % compute ngram vector for the testing window; the encoderMode determines the type of encoding algorithm   
        %global encoderMode;
		%switch encoderMode
            %case '4N_features'
		        %sigHV = computeNgram_4N_features (testSet (i : i+N-1,:), iM, D, N, percision);
                sigHV = computeNgram_4N_featuresSparse(testSet (i : i+N-1,:), iM, D, N, percision);
            %case 'N_feature_perm'
		        %sigHV = computeNgram_1N_feature_perm (testSet (i : i+N-1,:), iM, D, N, percision);
        %end

        % find the label of testing window
        maxAngle = -1;
        predicLabel = -1;
		for label = 1:1:max(labelTestSet)
			%angle = cosAngle(AM (label), sigHV);
            angle = sum(AM (label).* sigHV);
			if (angle > maxAngle)
				maxAngle = angle;
				predicLabel = label;
			end
        end
        %numpat(1,actualLabel)=numpat(1,actualLabel)+1;
		if predicLabel == actualLabel
			correct = correct + 1;
            %numpatpred(1,predicLabel)=numpatpred(1,predicLabel)+1;
            
        elseif labelTestSet (i) ~= labelTestSet(i+N-1)
			tranzError = tranzError + 1;
		end
    end
    
%     disp(numpat)
%     disp(numpatpred)
  
    accuracy = correct / numTests;
	accExcTrnz = (correct + tranzError) / numTests;
end

%% S



function [TR_data, TR_label] = featurizeData (labels, data,  N)
%
% DESCRIPTION   : test accuracy based on input testing data
%
% INPUTS:
%   labelTestSet: testing labels
%   testSet     : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy for all situations
%   accExcTrnz  : classification accuracy excluding the transitions between gestutes
%
            
            TR_data=[];
            TR_label=[];
            %for label = 1:1:max(labels)
				
			%	numPat (label) = 0;
			%end
    		i = 1;
            k=1;
			while i < length(labels)-N+1
				% if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
				if labels(i) == labels(i+N-1)		
                    
                    %numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
                    data_c=data(i : i+N-1,:);
                    TR_data(k,:)=data_c(:);
                    
                    TR_label(k,1)=labels(i+N-1);
                    k=k+1;

					i = i + 1;
				else
					% jump to the next window that has a 'stable' gesture where all labels are equal
					i = i + N - 1; %DK. double check
				end
            end
            %disp(k)
	%for label = 1:1:max(labels)
	%	fprintf ('Class= %d \t mean= %0.3f \t n_pat=%d \t created \n', label, numPat(label), numPat(label));
    %end
end


function [TS_data, TS_label] = featurizeTestingData (labelTestSet, testSet,  N)
%
% DESCRIPTION   : test accuracy based on input testing data
%
% INPUTS:
%   labelTestSet: testing labels
%   testSet     : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy for all situations
%   accExcTrnz  : classification accuracy excluding the transitions between gestutes
%
            
            TS_data=[];
            TS_label=[];
            numTests = 0;   
            k=1;

            	for i = 1:N:length(testSet)-N+1
                numTests = numTests + 1;
                % set actual label as the majority of labels in the window!
                actualLabel = mode(labelTestSet (i : i+N-1));
                data_c= testSet (i : i+N-1,:);
                TS_data(k,:)=data_c(:);  
                TS_label(k,1)=actualLabel;
                 k=k+1;
                % find the label of testing window
                
                end
           
            %for label = 1:1:max(labels)
				
			%	numPat (label) = 0;
			%end
            %disp(k)
	%for label = 1:1:max(labels)
	%	fprintf ('Class= %d \t mean= %0.3f \t n_pat=%d \t created \n', label, numPat(label), numPat(label));
    %end
end



function [iM] = initItemMemoriesHoloGN (D, MAXL, N,seed,density)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
% seed          : seed to set random generator in a prespecified position. This will randomize results due to initial vectors
% density       : set initial density of 1s for vectors.
% OUTPUTS:
%   iM          : item memory
 
	iM = containers.Map ('KeyType','char','ValueType','any');
    rng ('default');
    rng (seed);
       
    % Here we form N*4 features, hence we require N*4 orthogonal initial vectors        
    for j = 1:1:N*4
        currentHV =double(rand(1,D)>=(1-density));
                
        % Iterate over dicrete levels and 'continuously' generate vectors for related values	
        for i = 0:1:MAXL
            %disp(i)
            key = strcat(int2str(i),'x',int2str(j));
            iM(key) = currentHV;
            currentHV=circshift(currentHV, [0 1]);
        end
    end
end


%% Inititialize memory with nonlinear characteristic
function [iM] = initItemMemoriesNonLinear (D, MAXL, N,seed,density)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
% seed          : seed to set random generator in a prespecified position. This will randomize results due to initial vectors
% density       : set initial density of 1s for vectors.
% OUTPUTS:
%   iM          : item memory
 
	iM = containers.Map ('KeyType','char','ValueType','any');
    rng ('default');
    rng (seed);
    
    % Here we form N*4 features, hence we require N*4 orthogonal initial vectors        
    for j = 1:1:N*4
	    currentHV = genRandomHV (D,density);
	    %randomIndex = randperm (D);
        
        
        % Iterate over dicrete levels and 'continuously' generate vectors for related values	
        for i = 0:1:MAXL+40
            %disp(i)
            key = strcat(int2str(i),'x',int2str(j));
            iM(key) = currentHV;
       
            density_act=sum(currentHV==1); %actual number of ones
            pos1=find(currentHV==1);
            pos0=find(currentHV~=1);
            randomIndex1 = randperm (density_act);
            randomIndex0 = randperm (D-density_act);    
            k= floor((density_act-D*((density_act/D)^2))/MAXL);
            

                    
             
            
            global codingMode;
	        switch codingMode
                case 'dense_bipolar'
		            %currentHV (randomIndex(startInx : endInx)) = currentHV (randomIndex(startInx : endInx)) * -1;
                    %currentHV (pos1(randomIndex1(startInx : endInx))) = -1;
                    %currentHV (pos0(randomIndex0(startInx : endInx))) = 1;
                    currentHV (pos1(randomIndex1(1 : k))) = -1;
                    currentHV (pos0(randomIndex0(1 : k))) = 1;
                    
                case 'dense_binary'
		            %currentHV (randomIndex(startInx : endInx)) = not (currentHV (randomIndex(startInx : endInx)));
                    %currentHV (pos1(randomIndex1(startInx : endInx))) = 0;
                    %currentHV (pos0(randomIndex0(startInx : endInx))) = 1;
                    currentHV (pos1(randomIndex1(1 : k))) = 0;
                    currentHV (pos0(randomIndex0(1 : k))) = 1;                   
                    
            end
            
        end
    end
end

%% Inititialize memory with linear characteristic
function [iM] = initItemMemoriesLinear (D, MAXL, N,seed,density)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
% seed          : seed to set random generator in a prespecified position. This will randomize results due to initial vectors
% density       : set initial density of 1s for vectors.
% OUTPUTS:
%   iM          : item memory
 
	iM = containers.Map ('KeyType','char','ValueType','any');
    rng ('default');
    rng (seed);
    
    % Here we form N*4 features, hence we require N*4 orthogonal initial vectors        
    for j = 1:1:N*4
	    currentHV = genRandomHV (D,density);
        currentHV(2,:)=genRandomHV (D,density);
        
        

        
        % Iterate over dicrete levels and 'continuously' generate vectors for related values	
        for i = 0:1:MAXL
            %disp(i)
            key = strcat(int2str(i),'x',int2str(j));
            
            ind=round((1-(i/MAXL))*D);
            currentHV_lin=[currentHV(1,1:ind) ,currentHV(2,ind+1:D) ];
            iM(key) = currentHV_lin;                  
        end
    end
end

%% Inititialize projection matrix
function [iM] = initItemMemoriesProjection (D, MAXL, N,seed,density)
%
% DESCRIPTION   : initialize the item Memory (here continous) 
%
% INPUTS:
%   D           : Dimension of vectors
%   MAXL        : Maximum amplitude of EMG signal
%   N           : size of n-gram, i.e., window size
% seed          : seed to set random generator in a prespecified position. This will randomize results due to initial vectors
% density       : set initial density of 1s for vectors.
% OUTPUTS:
%   iM          : item memory
    rng ('default');
    rng (seed);
	
    iM=2*(round(rand(N*4,D))-0.5);
    %iM=randn(N*4,D);
   
end



%% 
function [numPat, AM_b] = hdctrainProjection (labels, data, iM, D, N, percision, cuttingAngle) 
%
% DESCRIPTION   : train an associative memory based on input training data
%
% INPUTS:
%   labels      : training labels
%   data        : EMG training data
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%   cuttingAngle: threshold angle for not including a vector into SUM vector
%
% OUTPUTS:
%   AM          : Trained associative memory
%   numPat      : Number of stored patterns for each class of AM
%
    %global Projected
    %Projected=[];
    
    %d_ng=0;% counter density of ngrams. will be used in majorutyCDT
    %c_ng=0;%counter for the number of ngrams. will be used in majoruty sum
    
    % initialize an empty AM
    AM = containers.Map ('KeyType','double','ValueType','any');
    AM_b = containers.Map ('KeyType','double','ValueType','any'); %binarized AM
	%fprintf ('Total traning samples size = %d\n', length(labels));
    for label = 1:1:max(labels)
    	AM (label) = zeros (1,D);
        AM_b(label) = zeros (1,D);
	    numPat (label) = 0;
    end

    % walk through the input data 
    i = 1;
    while i < length(labels)-N+1
        % if labels for the entire of n-gram (or window) are the same; it excludes training samples that are in transitions between 2 neighbor gestures 
        if labels(i) == labels(i+N-1)		
            % compute n-gram; the encoderMode determines the type of encoding algorithm
           % global encoderMode;
	        %switch encoderMode
                %case '4N_features'
                    
                    feature=data(i : i+N-1,:);
                    feature=feature(:)';
                    
                    feature_l=length(feature);
                    
                    
			        %ngram = computeNgram_4N_featuresSparse (data(i : i+N-1,:), iM, D, N, percision);
                    ngram =feature*iM(1:feature_l ,:) ;
                    %Projected(end+1,:)=ngram;
                    
               % case 'N_feature_perm'
			        %ngram = computeNgram_1N_feature_perm (data(i : i+N-1,:), iM, D, N, percision);
            %end

			% check whether the produced n-gram is already in AM
            %angle = cosAngle(ngram, AM (labels(i+N-1)));
            %angle = sum(AM_b (labels(i+N-1)).* ngram)/sum(ngram);

			%if angle < cuttingAngle | isnan(angle)
				AM (labels(i+N-1)) = AM (labels(i+N-1)) + ngram;
                numPat (labels(i+N-1)) = numPat (labels(i+N-1)) + 1;
                %d_ng=d_ng+sum(ngram);
                %c_ng=c_ng+1;
                %d_mean=round(d_ng/c_ng);
                %AM_b (labels(i+N-1)) =MajorityCDT( AM (labels(i+N-1)),d_mean);
                	            
			%end
            % walk through data with stride 1
            i = i + 1;
		else
            % jump to the next window that has a 'stable' gesture where all labels are equal
            i = i + N - 1;
        end
    end
        %d_mean=round(d_ng/c_ng);
        
    for label = 1:1:max(labels)
            AM_b (label) =AM (label);
            %AM_b (label) =MajorityCDT( AM (label),round(0.1*D)); %set a certain Percentage of the initial dimensionalty
            %disp(round(0.9*d_mean));
		%fprintf ('Class= %d \t mean= %.0f \t n_pat=%d \t created \n', label, mean(AM(label)), numPat(label));
    end
end


%% 
function [accExcTrnz, accuracy] = hdcpredictProjection (labelTestSet, testSet, AM, iM, D, N, percision)
%
% DESCRIPTION   : test accuracy based on input testing data
%
% INPUTS:
%   labelTestSet: testing labels
%   testSet     : EMG test data
%   AM          : Trained associative memory
%   iM          : item memory
%   D           : Dimension of vectors
%   N           : size of n-gram, i.e., window size
%   percision   : percision used in quantization of input EMG signals
%
% OUTPUTS:
%   accuracy    : classification accuracy for all situations
%   accExcTrnz  : classification accuracy excluding the transitions between gestutes
%
	correct = 0;        % count the number of correct classifications
    numTests = 0;       % count the total number of tests
	tranzError = 0;     % count the number of tests that went wrong but were between gesture transitions 
	%STAT=zeros(1,5);    %DK. To test simple hypotesis
    
    global testngrams;
    testngrams=[];
    k=1;
    
    % Go over the test data with stride of N
    for i = 1:N:length(testSet)-N+1
		numTests = numTests + 1;
		% set actual label as the majority of labels in the window!
        actualLabel = mode(labelTestSet (i : i+N-1));
        
        % compute ngram vector for the testing window; the encoderMode determines the type of encoding algorithm   
                    feature=testSet (i : i+N-1,:);
                    feature=feature(:)';
                    feature_l=length(feature);  
                    sigHV =feature*iM(1:feature_l ,:); %create a projection
                
                maxAngle = -1;
                predicLabel = -1;
                %maxAngle = Inf;
                %predicLabel = Inf;                
                
                for label = 1:1:max(labelTestSet)
                    
                    %correlation=corrcoef(AM (label),sigHV);
                    %angle=correlation(1,2);
                    angle = cosAngle(AM (label), sigHV);
                    testngrams(k,label)=angle;
                    
                    if (angle > maxAngle)
                        maxAngle = angle;
                        predicLabel = label;
                    end
                    
%                     angle =  sum( (AM (label)- sigHV).^2)^0.5;
%                     if (angle < maxAngle)
%                         maxAngle = angle;
%                         predicLabel < label;
%                     end
                    
                    
                    
                    
                end
                k=k+1;
        
                if predicLabel == actualLabel
                    correct = correct + 1;
                    
                elseif labelTestSet (i) ~= labelTestSet(i+N-1)
                    tranzError = tranzError + 1;
                end
        
        %STAT(1,actualLabel)=STAT(1,actualLabel)+1;
    end
    %disp(STAT)
    accuracy = correct / numTests;
	accExcTrnz = (correct + tranzError) / numTests;
end


%% 
    function [TR_data_c, TR_label_c] = make_centorid (TR_data, TR_label)
    
    TR_data_c=zeros(max(TR_label),length( TR_data(1,:)));
    TR_label_c=zeros(max(TR_label),1);
    
    for i=1:max(TR_label)
        TR_label_c(i,1)=i;
        TR_data_c(i,:)=  mean(TR_data((TR_label==i),:));        
    end
    
    
    end 
      

