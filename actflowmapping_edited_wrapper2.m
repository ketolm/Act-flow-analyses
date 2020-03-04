%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrapper for actflowmapping2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This version is for the     %
% second set of VOIs and such %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all;

load('mcvsaSubjects2.mat');
load('mcvsaContrastMatrix2.mat')

load('mcvsmSubjects2.mat');
load('mcvsmContrastMatrix2.mat')

load('mrestSubjects.mat');
load('mrestBetasMatrix.mat')

components = 8;

%% For a single task (mcvsa or mcvsm) and mrest
% change between mcvsa and mcvsm

match = 0;
% Checks mcvsaSubjects2 and mrestSubjects and finds where the mats meet up
for curTaskNum = 1:length(mcvsmSubjects2)   
    curTaskSub = mcvsmSubjects2(curTaskNum);
    for curRestNum = 1:length(mrestSubjects)
        curRestSub = mrestSubjects(curRestNum);
        if curTaskSub == curRestSub
            match = match + 1;
            subMatcher(match,:) = [curTaskNum curRestNum];   
            %break
        end
    end 
end
disp(subMatcher)


mcvsmCurated = mcvsmContrastMatrix2(:,:,subMatcher(:,1));
mrestCurated = mrestBetasMatrix(:,:,1,subMatcher(:,2));

[r_overall, p_overall, t_overall, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean, principalValues, indices, principalValError, pr_overall, pp_overall, pt_overall, pr_bytask, pp_bytask, principalPredMatrix, pr_bysubj, pr_avgfirst_bytask, pr_avgfirst_mean] = actflowmapping_edited(mcvsmCurated, mrestCurated, components);

%% Split PD and Control Single Task

for curSubNum = 1:length(subMatcher)
    
    fileID = fopen(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mcvsmSubjects2(subMatcher(curSubNum))), '/session1/0_group'), 'r');
    if fileID == -1
        disp([int2str(mcvsmSubjects2(subMatcher(curSubNum))) 'doesnt have 0_group']);
    else
        diseaseStatus = fscanf(fileID, '%s');
        if strcmp(diseaseStatus,'PD')
            subDisease(curSubNum) = 1;
        else
            subDisease(curSubNum) = 0;
        end
    end
end

PD = 0;
if PD
    % for PD
    mcvsmPDSubjs = subMatcher(:,1) & subDisease';
    mcvsmPDSubjs = subMatcher(mcvsmPDSubjs,1);
    
    mrestPDsubjs = subMatcher(:,2) & subDisease';
    mrestPDsubjs = subMatcher(mrestPDsubjs,2);
    
    curatedTask  = mcvsmContrastMatrix2(:,:,mcvsmPDSubjs);
    mrestCurated = mrestBetasMatrix(:,:,1,mrestPDsubjs);
    
else
    % for controls
    mcvsmConSubjs = subMatcher(:,1) & (subDisease ~= 1)';
    mcvsmConSubjs = subMatcher(mcvsmConSubjs,1);
    
    mrestConSubjs = subMatcher(:,2) & (subDisease ~= 1)';
    mrestConSubjs = subMatcher(mrestConSubjs,2);
    
    curatedTask  =  mcvsmContrastMatrix2(:,:,mcvsmConSubjs);
    mrestCurated = mrestBetasMatrix(:,:,1,mrestConSubjs);
end

[r_overall, p_overall, t_overall, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean, principalValues, indices, principalValError, pr_overall, pp_overall, pt_overall, pr_bytask, pp_bytask, principalPredMatrix, pr_bysubj, pr_avgfirst_bytask, pr_avgfirst_mean] = actflowmapping_edited(curatedTask, mrestCurated,components);

%% For mcvsa + mcvsm and mrest

match = 0;
subMatcher = [0 0 0];
% Checks mcvsaSubjects2 and mrestSubjects and finds where the mats meet up
for curMCVSANum = 1:length(mcvsaSubjects2)   
    curMCVSASub = mcvsaSubjects2(curMCVSANum);
    
    for curMCVSMNum = 1:length(mcvsmSubjects2)
        curMCVSMSub = mcvsmSubjects2(curMCVSMNum);
        
        for curRestNum = 1:length(mrestSubjects)
            curRestSub = mrestSubjects(curRestNum);
            
            if curRestSub == curMCVSASub && curRestSub == curMCVSMSub

                match = match + 1;
                subMatcher(match,:) = [curMCVSANum curMCVSMNum curRestNum];   
                break
            end
        end 
    end
end
%disp(subMatcher)

curatedTask = cat(3, mcvsaContrastMatrix2(:,:,subMatcher(:,1)), mcvsmContrastMatrix2(:,1,subMatcher(:,2)));
mrestCurated = mrestBetasMatrix(:,:,1,subMatcher(:,3));

[r_overall, p_overall, t_overall, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean, principalValues, indices, principalValError, pr_overall, pp_overall, pt_overall, pr_bytask, pp_bytask, principalPredMatrix, pr_bysubj, pr_avgfirst_bytask, pr_avgfirst_mean] = actflowmapping_edited(curatedTask, mrestCurated, components);

%% Split PD and Control Combined

for curSubNum = 1:length(subMatcher)
    
    fileID = fopen(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mcvsaSubjects2(subMatcher(curSubNum))), '/session1/0_group'), 'r');
    if fileID == -1
        disp([int2str(mcvsaSubjects2(subMatcher(curSubNum))) 'doesnt have 0_group']);
    else
        diseaseStatus = fscanf(fileID, '%s');
        if strcmp(diseaseStatus,'PD')
            subDisease(curSubNum) = 1;
        else
            subDisease(curSubNum) = 0;
        end
    end
end

PD = 0;
if PD
    % for PD
    mcvsaPDSubjs = subMatcher(:,1) & subDisease';
    mcvsaPDSubjs = subMatcher(mcvsaPDSubjs,1);
    
    mcvsmPDsubjs = subMatcher(:,2) & subDisease';
    mcvsmPDsubjs = subMatcher(mcvsmPDsubjs,2);
    
    mrestPDsubjs = subMatcher(:,3) & subDisease';
    mrestPDsubjs = subMatcher(mrestPDsubjs,3);
    
    curatedTask = cat(3, mcvsaContrastMatrix2(:,:,mcvsaPDSubjs), mcvsmContrastMatrix2(:,1,mcvsmPDsubjs));
    mrestCurated = mrestBetasMatrix(:,:,1,mrestPDsubjs);
    
else
    % for controls
    mcvsaConSubjs = subMatcher(:,1) & (subDisease ~= 1)';
    mcvsaConSubjs = subMatcher(mcvsaConSubjs,1);
    
    mcvsmConSubjs = subMatcher(:,2) & (subDisease ~= 1)';
    mcvsmConSubjs = subMatcher(mcvsmConSubjs,2);
    
    mrestConSubjs = subMatcher(:,3) & (subDisease ~= 1)';
    mrestConSubjs = subMatcher(mrestConSubjs,3);
    
    curatedTask = cat(3, mcvsaContrastMatrix2(:,:,mcvsaConSubjs), mcvsmContrastMatrix2(:,1,mcvsmConSubjs));
    mrestCurated = mrestBetasMatrix(:,:,1,mrestConSubjs);
end

all_r_overall   = zeros(263,1);
all_r_avgfirst  = zeros(263,1);
all_pr_overall  = zeros(263,1);
all_pr_avgfirst = zeros(263,1);
    
for curNumComp = 1:263
    components = curNumComp;
    [r_overall, p_overall, t_overall, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean, principalValues, indices, principalValError, pr_overall, pp_overall, pt_overall, pr_bytask, pp_bytask, principalPredMatrix, pr_bysubj, pr_avgfirst_bytask, pr_avgfirst_mean] = actflowmapping_edited(curatedTask, mrestCurated, components);
    all_r_overall(curNumComp) = r_overall;
    all_r_avgfirst(curNumComp) = r_avgfirst_mean;
    all_pr_overall(curNumComp)  = pr_overall;
    all_pr_avgfirst(curNumComp) = pr_avgfirst_mean;
end

 plot([1:262],all_pr_avgfirst)
 xlim([1 264])
 ylim([.60 .95])

%% Extra shit looking at how numbers of components relates to r and p-values

match = 0;
subMatcher = [0 0 0];
% Checks mcvsaSubjects2 and mrestSubjects and finds where the mats meet up
for curMCVSANum = 1:length(mcvsaSubjects2)   
    curMCVSASub = mcvsaSubjects2(curMCVSANum);
    
    for curMCVSMNum = 1:length(mcvsmSubjects2)
        curMCVSMSub = mcvsmSubjects2(curMCVSMNum);
        
        for curRestNum = 1:length(mrestSubjects)
            curRestSub = mrestSubjects(curRestNum);
            
            if curRestSub == curMCVSASub && curRestSub == curMCVSMSub

                match = match + 1;
                subMatcher(match,:) = [curMCVSANum curMCVSMNum curRestNum];   
                break
            end
        end 
    end
end
%disp(subMatcher)

curatedTask = cat(3, mcvsaContrastMatrix2(:,:,subMatcher(:,1)), mcvsmContrastMatrix2(:,1,subMatcher(:,2)));
mrestCurated = mrestBetasMatrix(:,:,1,subMatcher(:,3));

all_r_overall   = zeros(263,1);
all_r_avgfirst  = zeros(263,1);
all_pr_overall  = zeros(263,1);
all_pr_avgfirst = zeros(263,1);
    
for curNumComp = 1:263
    components = curNumComp;
    [r_overall, p_overall, t_overall, r_bytask, p_bytask, taskActualMatrix, taskPredMatrix, r_bysubj, r_avgfirst_bytask, r_avgfirst_mean, principalValues, indices, principalValError, pr_overall, pp_overall, pt_overall, pr_bytask, pp_bytask, principalPredMatrix, pr_bysubj, pr_avgfirst_bytask, pr_avgfirst_mean] = actflowmapping_edited(curatedTask, mrestCurated, components);
    all_r_overall(curNumComp) = r_overall;
    all_r_avgfirst(curNumComp) = r_avgfirst_mean;
    all_pr_overall(curNumComp)  = pr_overall;
    all_pr_avgfirst(curNumComp) = pr_avgfirst_mean;
end

%% Looking at indices of the most important parts

principalIndices = zeros(264, 264, size(indices,4));


for i = 1:size(indices,1)  
    for j = 1:size(indices,4)
        for k = indices(i,:,:,j)
            principalIndices(i,k,j) = principalIndices(i,k,j) + 1;          
        end
    end
end

figure;
imshow(principalIndices(:,:,3))

figure;
imshow(mean(principalIndices,3))

meanP = mean(principalIndices,3);
sumP  = sum(meanP);
hist(sumP)
bar(1:264,sumP)
xlim([0 270])

bar(1:264,mean(meanP))

%%

[pRows, pCols] = find(sumP > 10.4);
pCoors = [pRows pCols];