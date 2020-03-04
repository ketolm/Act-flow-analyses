% Baddeley's Delta Metric
clc; clear all; close all;

% load the matrices
load('controlConnectivity')
load('PDConnectivity')


%% Base Level Metric

% binary threshold 
bThresh = .9;

% min distance threshold
cThresh = 7;

numCon = size(controlConnectivity,3);
numPD  = size(PDConnectivity,3);
controlBW  = nan(264,264,numCon);
controlEuc = nan(264,264,numCon);
PDBW  = nan(264,264,numCon);
PDEuc = nan(264,264,numCon);


% convert to binary images based on threshold
% create point-to-set distance image
for curNum = 1:numCon      
    controlBW(:,:,curNum) = squeeze(controlConnectivity(:,:,curNum)) > bThresh;
    allConDist = bwdist(controlBW(:,:,curNum));
    controlEuc(:,:,curNum) = allConDist .* (allConDist > cThresh);
    
    PDBW(:,:,curNum) = squeeze(PDConnectivity(:,:,curNum)) > bThresh;
    allPDDist = bwdist(PDBW(:,:,curNum));
    PDEuc(:,:,curNum) = allPDDist .* (allPDDist > cThresh);
end


% calculate absolute difference across all permutations
combos = nchoosek(1:numCon,2);
for n = 1:length(combos)
    curDiff = imabsdiff(controlEuc(:,:,combos(n,1)), PDEuc(:,:,combos(n,2)));
    diff(:,:,n) = curDiff;
end
%imshow(diff(:,:,4))

meanDiff = mean(diff,3);
imshow(meanDiff)

normDiff = meanDiff - min(meanDiff(:));
normDiff = normDiff ./ max(normDiff(:));
%imshow(normDiff)
%imshow(im2bw(normDiff))



%% Permutation Tests

perm = true;
permutationNum = 1000;

if perm
    maxDiff = nan(1,permutationNum);
    
    for curPerm = 1:permutationNum

        scramble  = cat(3, controlConnectivity, PDConnectivity(:,:,1:numCon));
        
        randOrder = randperm(size(scramble,3));
        scramble  = scramble(:,:,randOrder);
        
        randAssign = randperm(numCon * 2);
        rGroup1 = scramble(:,:,randAssign(1:numCon));
        rGroup2 = scramble(:,:,randAssign(numCon+1:end));
        
        rGroup1BW  = nan(264,264,numCon);
        rGroup1Euc = nan(264,264,numCon);
        rGroup2BW  = nan(264,264,numCon);
        rGroup2Euc = nan(264,264,numCon);
        
        
        for curNum = 1:numCon      
            rGroup1BW(:,:,curNum) = squeeze(rGroup1(:,:,curNum)) > bThresh;
            allGroup1Dist = bwdist(rGroup1BW(:,:,curNum));
            rGroup1Euc(:,:,curNum) = allGroup1Dist .* (allGroup1Dist > cThresh);
            
            rGroup2BW(:,:,curNum) = squeeze(rGroup2(:,:,curNum)) > bThresh;
            allGroup2Dist = bwdist(rGroup2BW(:,:,curNum));
            rGroup2Euc(:,:,curNum) = allGroup2Dist .* (allGroup2Dist > cThresh);
        end
        
        combos = nchoosek(1:numCon,2);
        for n = 1:length(combos)
            curDiff = imabsdiff(rGroup1Euc(:,:,combos(n,1)), rGroup2Euc(:,:,combos(n,2)));
            rDiff(:,:,n) = curDiff;
        end
        
        rMeanDiff = mean(rDiff,3);
        
        maxDiff(curPerm) = max(max(rMeanDiff));
        
    end
end

%% Visualize 
close all;

figure;
hist(maxDiff)
sortedMaxDiff = sort(sort(maxDiff));
sigThresh = sortedMaxDiff(permutationNum * .95);

figure;
sigMeanDiff = meanDiff > sigThresh;
imshow(sigMeanDiff)

%%
[xs, ys] = find(sigMeanDiff);
coors = [xs ys];

%sigValues = meanDiff .* sigMeanDiff;

for curCoor = 1:size(coors)
    
    sigValues(curCoor) = meanDiff(coors(curCoor,1), coors(curCoor,2));
    
end
