%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences in beta connectivity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Didn't find shit %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all;

load('mrestSubjects.mat');

for curSubNum = 1:length(mrestSubjects)
    
    fileID = fopen(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mrestSubjects(curSubNum)), '/session1/0_group'), 'r');
    if fileID == -1
        disp([int2str(mrestSubjects(curSubNum)) 'doesnt have 0_group']);
    else
        diseaseStatus = fscanf(fileID, '%s');
        if strcmp(diseaseStatus,'PD')
            subDisease(curSubNum) = 1;
        else
            subDisease(curSubNum) = 0;
        end
    end
end


numPD = 0;
numControl = 0;

for curNum = 1:length(mrestSubjects)
    
    if subDisease(curNum) == 1
        numPD = numPD + 1;
        load(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mrestSubjects(curNum)), '/session1/mrest_results/activity_flow/netMat_betas.mat'))
        PDConnectivity(:,:,numPD) = netMat_betas;
    else
        numControl = numControl + 1;
        load(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mrestSubjects(curNum)), '/session1/mrest_results/activity_flow/netMat_betas.mat'))
        controlConnectivity(:,:,numControl) = netMat_betas;
    end
end

%% MSE Within Groups

combos = nchoosek(1:numControl,2);
for n = 1:length(combos)
    % homegrown MSE
    PDArrayMSE(n)       = sum(sum((PDConnectivity(:,:,combos(n,1)) - PDConnectivity(:,:,combos(n,2))).^2)) / (264^2 - 264);
    controlArrayMSE(n)  = sum(sum((controlConnectivity(:,:,combos(n,1)) - controlConnectivity(:,:,combos(n,2))).^2)) / (264^2 - 264);
    
    % MATLAB MSE
    [PDPSNR(n), PDMSE(n), PDMAXERR(n), PDL2RAT(n)]                     = measerr(PDConnectivity(:,:,combos(n,1)) , PDConnectivity(:,:,combos(n,2)));
    [controlPSNR(n), controlMSE(n), controlMAXERR(n), controlL2RAT(n)] = measerr(controlConnectivity(:,:,combos(n,1)) , controlConnectivity(:,:,combos(n,2)));
end

mPDMSE    = mean(PDArrayMSE);
mPDMSE2   = mean(PDMSE);
mPDPSNR   = mean(PDPSNR);
mPDMAXERR = mean(PDMAXERR);
mPDL2RAT  = mean(PDL2RAT);

mControlMSE    = mean(controlArrayMSE);
mControlMSE2   = mean(controlMSE);
mControlPSNR   = mean(controlPSNR);
mControlMAXERR = mean(controlMAXERR);
mControlL2RAT  = mean(controlL2RAT);



%% MSE Between Groups

MSE  = 0;
PSNR = 0;
MSE2 = 0;
MAXERR = 0;
L2RAT = 0;

combos = nchoosek(1:numControl,2);
for n = 1:length(combos)
    MSE(n)  = sum(sum((PDConnectivity(:,:,combos(n,1)) - controlConnectivity(:,:,combos(n,2))).^2)) / (264^2 - 264);
    [PSNR(n), MSE2(n), MAXERR(n), L2RAT(n)] = measerr(PDConnectivity(:,:,combos(n,1)) , controlConnectivity(:,:,combos(n,2)));
end

mMSE    = mean(MSE);
mMSE2   = mean(MSE2);
mPSNR   = mean(PSNR);
mMAXERR = mean(MAXERR);
mL2RAT  = mean(L2RAT);


%% Permutation testing

permutation = 0;

if permutation

    groupMembership = shuffles([zeros(1,numControl) ones(1,numControl)]);

    n = 10000;

    mRandMSE    = ones(1,n)*-1;
    mRandMSE2   = ones(1,n)*-1;
    mRandPSNR   = ones(1,n)*-1;
    mRandMAXERR = ones(1,n)*-1;
    mRandL2RAT  = ones(1,n)*-1;

    for iteration = 1:n

        groupMembership = shuffles(groupMembership);

        numPD = 0;
        numControl = 0;

        for curSub = 1:length(groupMembership)

            if groupMembership(curSub) == 1
                numPD = numPD + 1;
                load(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mrestSubjects(curSub)), '/session1/mrest_results/activity_flow/netMat_betas.mat'))
                randPDConnectivity(:,:,numPD) = netMat_betas;
            else
                numControl = numControl + 1;
                load(strcat('/mnt/praxic/pdnetworks2/subjects/', int2str(mrestSubjects(curSub)), '/session1/mrest_results/activity_flow/netMat_betas.mat'))
                randControlConnectivity(:,:,numControl) = netMat_betas;
            end
        end


        for curCom = 1:length(combos)
            randMSE(curCom)  = sum(sum((randPDConnectivity(:,:,combos(curCom,1)) - randControlConnectivity(:,:,combos(curCom,2))).^2)) / (264^2 - 264);
            [randPSNR(curCom), randMSE2(curCom), randMAXERR(curCom), randL2RAT(curCom)] = measerr(randPDConnectivity(:,:,combos(curCom,1)), randControlConnectivity(:,:,combos(curCom,2)));
        end

        mRandMSE(iteration)    = mean(randMSE);
        mRandMSE2(iteration)   = mean(randMSE2);
        mRandPSNR(iteration)   = mean(randPSNR);
        mRandMAXERR(iteration) = mean(randMAXERR);
        mRandL2RAT(iteration)  = mean(randL2RAT);

    end
end



