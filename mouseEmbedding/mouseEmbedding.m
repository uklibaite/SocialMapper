%% make mouse embedding

sDANNCE_T = readtable('sDANNCE_file_info.xlsx');
cohort = table2cell(sDANNCE_T(:,1));
[groups, fid, groupidall] = unique(cohort);
mouseidx = find(groupidall==9);

allFiles = table2cell(sDANNCE_T(:,2));
allmouseFiles = allFiles(mouseidx);

% to preserve order here for downstream processing - use these paths
% instead

load('mouseFileOrder.mat')


%% load all mouse data
% downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VKJHTD
% cd dataverse_files

allMOUSEL = cell(32,1);
for i = 1:32
    load(mouseOrderShort{i});
    % sdannce (struct with info)
    allMOUSEL{i} = sdannce.m1;
end

allMOUSES = cell(48,2);
for i = 1:48
    xi = i+32;
    load(mouseOrderShort{i});
    allMOUSES{i,1} = sdannce.m1;
    allMOUSES{i,2} = sdannce.m2;
end


%% look at data

% skeleton - load from skeletons/rat23.mat 
joints_idx = [1 2; 1 3; 2 3; ...
    1 4; 4 5; 5 6; 6 7; ...
    4 8; 8 9; 9 10; 10 11; ...
    4 12; 12 13; 13 14; 14 15; ...
    6 16; 16 17; 17 18; 18 19; ...
    6 20; 20 21; 21 22; 22 23];

chead = [1 .6 .2]; % orange
% cspine = [198 252 3]./256; % yellow
% cspine = [252 215 3]./256; % yellow
cspine = [.2 .635 .172]; % green
cLF = [0 0 1]; % blue
cRF = [1 0 0]; % red
cLH = [0 1 1]; % cyan
cRH = [1 0 1]; % magenta
scM = [chead; chead; chead; cspine; cspine; cspine; cspine; cLF; cLF; ...
    cLF; cLF; cRF; cRF; cRF; cRF;...
    cLH; cLH; cLH; cLH; cRH; cRH; cRH; cRH];
sc = scM;
skeleton.color = scM;%zeros(23,3); 
skeleton.joints_idx = joints_idx;


% single animal:

close all
fIndic= figure('Name','Lone Test');
h = cell(1);
h{1} = Keypoint3DAnimator(sdannce.m1, skeleton ,'MarkerSize',15);
set(gca,'Color','w')
Animator.linkAll(h)
axis equal;
axis([-200 250 -150 350 -10 150]);
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);


% social:

psoc = cat(3,sdannce.m1, sdannce.m2);
skelSOC.joints_idx = [skeleton.joints_idx; skeleton.joints_idx+23];
skelSOC.mcolor = [sc; sc];
skelSOC.color = zeros(46,3); skelSOC.color(1:23,:) = .5*ones(23,3);
skelSOC2 = skelSOC;
skelSOC2.color(1:23,:) = ones(23,3);
skelSOC2.color(24:end,:) = .7*ones(23,3);

skelSOC.mcolor = zeros(46,3);
skelSOC.color(1:23,:) = repmat([0 0 1], [23 1]); %ones(23,3);
skelSOC.color(24:end,:) = repmat([1 0 0], [23 1]); %.7*ones(23,3);


close all
fIndic= figure('Name','Social Test');
h = cell(1);
h{1} = Keypoint3DAnimator(psoc, skelSOC, 'Position', [0 0  1 1],'MarkerSize',15);
set(gca,'Color','w')
Animator.linkAll(h)
axis equal;
axis([-200 250 -150 350 -10 150]);
set(gcf,'Color','white');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);




%% LONE EMBEDDING

xIdx = 1:23; yIdx = 1:23;
[Xi Yi] = meshgrid(xIdx,yIdx);
Xi = Xi(:); Yi = Yi(:);
IDX = find(Xi~=Yi);
nx = length(xIdx);
firstBatch = true;
currentImage = 0;
batchSize = 30000;

% PCA
allMOUSE = [allMOUSEL; allMOUSES(:)];
lengtht = [lengthtL; lengthtS(:)];
for j = 1:length(allMOUSE)
    try
        fprintf(1,['Processing batch # ' num2str(j) '\n']);
        ma1 = allMOUSE{j};
        nn1 = size(ma1,1);
        p1Dist = zeros(nx^2,size(ma1,1));
        for i = 1:size(p1Dist,1)
            p1Dist(i,:) = returnDist3d(squeeze(ma1(:,:,Xi(i))),squeeze(ma1(:,:,Yi(i))));
        end
        p1Dsmooth = zeros(size(p1Dist));
        for i = 1:size(p1Dist,1)
            p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
        end
        p1Dist = p1Dsmooth(IDX,:)';
        
        scaleVal = lengtht(j)./90;
        p1Dist = p1Dist.*scaleVal;
        
        if firstBatch
            firstBatch = false;
            if size(p1Dist,1) < batchSize
                cBatchSize = size(p1Dist,1);
                X = p1Dist;
            else
                cBatchSize = batchSize;
                X = p1Dist;
            end
            currentImage = cBatchSize;
            mu = sum(X);
            C = cov(X).*cBatchSize + (mu'*mu)./ cBatchSize;
        else
            if size(p1Dist,1) < batchSize
                cBatchSize = size(p1Dist,1);
                X = p1Dist;
            else
                cBatchSize = batchSize;
                X = p1Dist(randperm(size(p1Dist,1),cBatchSize),:);
            end
            tempMu = sum(X);
            mu = mu + tempMu;
            C = C + cov(X).*cBatchSize + (tempMu'*tempMu)./cBatchSize;
            currentImage = currentImage + cBatchSize;
        end
    catch
    end
end

L = currentImage; mu = mu ./ L; C = C ./ L - mu'*mu;
fprintf(1,'Finding Principal Components\n');
[vecs,vals] = eig(C); vals = flipud(diag(vals)); vecs = fliplr(vecs);
mus = mu;
% save('vecsValsMouse0519.mat','C','L','mus','vals','vecs');

vecs15 = vecs(:,1:15);
minF = .5; maxF = 20; pcaModes = 20; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = 320;


mD = cell(size(allMOUSE)); mA = cell(size(allMOUSE));
for j = 1:length(allMOUSE)
    try
    j
    clear ma ma1 m1Dist nnData 
    ma1 = allMOUSE{j};

    nn1 = size(ma1,1);
    p1Dist = zeros(nx^2,size(ma1,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(ma1(:,:,Xi(i))),squeeze(ma1(:,:,Yi(i))));
    end

    p1Dist = p1Dist(IDX,:)';
    allz = squeeze(ma1(:,3,[19 23])); zz = allz(:);
    fz = prctile(zz,10);
    sj = returnDist3d(squeeze(ma1(:,:,1)),squeeze(ma1(:,:,7)));
    lz = prctile(sj,95);
    scaleVal = 90./lz;
    p1Dist = p1Dist.*scaleVal;
    
    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs15;
   
    [data,~] = findWavelets(proj,numModes,parameters);
    
    n = size(p1Dist,1);
    amps = sum(data,2);
    data2 = log(data);
    data2(data2<-5) = -5;
    jv = zeros(n,length(xIdx));
    for i = 1:length(xIdx)
        jv(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma1(:,:,xIdx(i)))).^2,2)),10)];
    end
    jv = jv.*scaleVal;
    jv(jv>=5) = 5;
    
    p1z = zeros(nx,nn1);
    for i = 1:nx
        p1z(i,:) = smooth(medfilt1(squeeze(ma1(:,3,xIdx(i))),3),3);
    end
    allz1 = squeeze(ma1(:,3,[19 23])); zz1 = allz1(:); fz1 = prctile(zz1,10);
    floorval = fz1;
    p1z = (p1z-floorval).*scaleVal; 
    
    nnData = [data2 .25*p1z' .5*jv];
    fprintf(1,'Running tSNE \n');
    tic
    yData = tsne(nnData(1:20:end,:));
    toc
    [signalData,signalAmps] = findTemplatesFromData(...
        nnData(1:20:end,:),yData,amps(1:20:end,:),numPerDataSet,parameters);
    mD{j} = signalData; mA{j} = signalAmps;
    catch 
    end
end
% save('loneSignalDataAmps_0519.mat','mA','mD');

allD = combineCells(mD,1); allA = combineCells(mA);
for i = 1:10
    tic
    Y{i} = tsne(allD);
    toc
end

yData = Y{7};
% save('train0519.mat','yData','allD')


%% Reembeddings lone
load('vecsValsMouse0519.mat','C','L','mus','vals','vecs');
load('train0519.mat','yData','allD');
trainingSetData = allD; trainingEmbeddingZ = yData;

vecs15 = vecs(:,1:15);
minF = .5; maxF = 20; pcaModes = 20; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.batchSize = 10000;

for j = 1:length(allMOUSE)
    clear ma ma1 m1Dist nnData
    ma1 = allMOUSE{j};

    nn1 = size(ma1,1);
    p1Dist = zeros(nx^2,size(ma1,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(ma1(:,:,Xi(i))),squeeze(ma1(:,:,Yi(i))));
    end

    p1Dist = p1Dist(IDX,:)';
    allz = squeeze(ma1(:,3,[19 23])); zz = allz(:);
    fz = prctile(zz,10);
    sj = returnDist3d(squeeze(ma1(:,:,1)),squeeze(ma1(:,:,7)));
    lz = prctile(sj,95);
    scaleVal = 90./lz;
    p1Dist = p1Dist.*scaleVal;

    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs15;

    [data,~] = findWavelets(proj,numModes,parameters);

    n = size(p1Dist,1);
    amps = sum(data,2);
    data2 = log(data);x
    data2(data2<-5) = -5;
    jv = zeros(n,length(xIdx));
    for i = 1:length(xIdx)
        jv(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma1(:,:,xIdx(i)))).^2,2)),10)];
    end
    jv = jv.*scaleVal;
    jv(jv>=5) = 5;

    p1z = zeros(nx,nn1);
    for i = 1:nx
        p1z(i,:) = smooth(medfilt1(squeeze(ma1(:,3,xIdx(i))),3),3);
    end
    allz1 = squeeze(ma1(:,3,[19 23])); zz1 = allz1(:); fz1 = prctile(zz1,10);
    floorval = fz1;
    p1z = (p1z-floorval).*scaleVal;

    nnData = [data2 .25*p1z' .5*jv];

    fprintf(1,'Finding Embeddings\n');
    [zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
        findTDistributedProjections_fmin(nnData,trainingSetData,...
        trainingEmbeddingZ,[],parameters);

    z = zValues; z(~inConvHull,:) = zGuesses(~inConvHull,:);
    filename = allMOUSE_files{i}; 
    % save(['mouse/RE_lone/RE_LONE_' num2str(j) '.mat'],'z','inConvHull','filename')
end




%%
useID = [1 4 5 6 7 8 11 12 15 16 19 20 23];
xIdx = [useID useID+23];
yIdx = [useID useID+23];

[Xi Yi] = meshgrid(xIdx,yIdx);
Xi = Xi(:); Yi = Yi(:);
IDX = find(Xi~=Yi);
nx = length(xIdx);
firstBatch = true;
currentImage = 0;
batchSize = 90000;

minF = 1/5; maxF = 5; pcaModes = 69; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.numPeriods = 25;

% PCA

for rec = 1:length(allMOUSES)
    try
        fprintf(1,['Processing batch # ' num2str(rec) '\n']);

        ma1 = allMOUSES{rec,1};
        ma2 = allMOUSES{rec,2};
        maS = cat(3,ma1,ma2);

        allz = squeeze(maS(:,3,[19 23 42 46])); zz = allz(:);
        fz = prctile(zz,10);

        nn1 = size(maS,1);
        p1Dist = zeros(nx^2,size(maS,1));
        for i = 1:size(p1Dist,1)
            p1Dist(i,:) = returnDist3d(squeeze(maS(:,:,Xi(i))),squeeze(maS(:,:,Yi(i))));
        end
        p1Dsmooth = zeros(size(p1Dist));
        for i = 1:size(p1Dist,1)
            p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
        end
        p1Dist = p1Dsmooth(IDX,:)';

        if firstBatch
            firstBatch = false;
            if size(p1Dist,1) < batchSize
                cBatchSize = size(p1Dist,1);
                X = p1Dist;
            else
                cBatchSize = batchSize;
                X = p1Dist;
            end
            currentImage = cBatchSize;
            mu = sum(X);
            C = cov(X).*cBatchSize + (mu'*mu)./ cBatchSize;
        else
            if size(p1Dist,1) < batchSize
                cBatchSize = size(p1Dist,1);
                X = p1Dist;
            else
                cBatchSize = batchSize;
                X = p1Dist(randperm(size(p1Dist,1),cBatchSize),:);
            end
            tempMu = sum(X);
            mu = mu + tempMu;
            C = C + cov(X).*cBatchSize + (tempMu'*tempMu)./cBatchSize;
            currentImage = currentImage + cBatchSize;
        end
    catch
    end
end
L = currentImage; mu = mu ./ L; C = C ./ L - mu'*mu;
fprintf(1,'Finding Principal Components\n');
[vecs,vals] = eig(C); vals = flipud(diag(vals)); vecs = fliplr(vecs);
mus = mu;
% save('vecsValsSOCMouse0520.mat','C','L','mus','vals','vecs');

vecs6 = vecs(:,1:6);
minF = .2; maxF = 5; pcaModes = 10; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = 400;


mDC = cell(size(allMOUSES)); mDCA = cell(size(allMOUSES));
for rec = 1:length(allMOUSES)
    rec
    ma1 = allMOUSES{rec,1};
    ma2 = allMOUSES{rec,2};
    maS = cat(3,ma1,ma2);
    allz = squeeze(maS(:,3,[19 23 42 46])); zz = allz(:);
    fz = prctile(zz,10);

    nn1 = size(maS,1);
    p1Dist = zeros(nx^2,size(maS,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(maS(:,:,Xi(i))),squeeze(maS(:,:,Yi(i))));
    end
    p1Dsmooth = zeros(size(p1Dist));
    for i = 1:size(p1Dist,1)
        p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
    end
    p1Dist = p1Dsmooth(IDX,:)';

    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs6;

    % velocity information
    p1v = zeros(nn1,3); p2v = zeros(nn1,3);
    for i = 1:23
        p1v(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma1(:,:,i))).^2,2)),10)];
        p2v(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma2(:,:,i))).^2,2)),10)];
    end

    p1v(p1v>5) = 5; p2v(p2v>5) = 5;
    n1 = ma1(:,:,1); s1 = ma1(:,:,4); c1 = ma1(:,:,5); t1 = ma1(:,:,7);
    n2 = ma2(:,:,1); s2 = ma2(:,:,4); c2 = ma2(:,:,5); t2 = ma2(:,:,7);
    distc = sqrt(sum((c1-c2).^2,2));

    % angle information
    angInfo = zeros(nn1,6);
    vec1 = n2-s1; vec2 = n1-s1; angInfo(:,1) = findAngleDiff(vec1,vec2);
    vec1 = c2-s1; vec2 = n1-s1; angInfo(:,2) = findAngleDiff(vec1,vec2);
    vec1 = t2-s1; vec2 = n1-s1; angInfo(:,3) = findAngleDiff(vec1,vec2);
    vec1 = n1-s2; vec2 = n2-s2; angInfo(:,4) = findAngleDiff(vec1,vec2);
    vec1 = c1-s2; vec2 = n2-s2; angInfo(:,5) = findAngleDiff(vec1,vec2);
    vec1 = t1-s2; vec2 = n2-s2; angInfo(:,6) = findAngleDiff(vec1,vec2);
    angInfo = pi - angInfo;

    d2n = sqrt(sum((n1-n2).^2,2)); % dist of nose1 to other rat
    d2t = sqrt(sum((n1-t2).^2,2));
    d2c = sqrt(sum((n1-c2).^2,2));
    d2anyt = [d2n d2c d2t];
    d2anyN = d2anyt./distc;

    d1n = sqrt(sum((n2-n1).^2,2)); % dist of nose2 to other rat
    d1t = sqrt(sum((n2-t1).^2,2));
    d1c = sqrt(sum((n2-c1).^2,2));
    d1anyt = [d1n d1c d1t];
    d1anyN = d1anyt./distc;

    PJ = [d1anyt d2anyt];
    PJN = [d1anyN d2anyN];
    [W,~] = findWavelets(proj,9,parameters);
    data2 = log(W); data2(data2<-3) = -3;

    % compute other social features and make final feature matrix
    rat1z = squeeze(ma1(:,3,[1 4 5 8 11 12 15]))-fz;
    rat2z = squeeze(ma2(:,3,[1 4 5 8 11 12 15]))-fz;
    nnData = [data2 .25*rat1z .25*rat2z 2*p1v(:,[1 11 12 15 19 23]) 2*p2v(:,[1 11 12 15 19 23]) 10*angInfo .05*PJ 2*PJN];
    yData = tsne(nnData(1:20:end,:));
    amps = sum(nnData',1);
    dataMat1 = nnData;

    % FIND GOOD TEMPLATES AND SAVE THEM
    [signalData,signalAmps] = findTemplatesFromData(...
        nnData(1:20:end,:),yData,amps(1:20:end,:),numPerDataSet,parameters);
    mDC{rec,1} = signalData; mDCA{rec,1} = signalAmps;

    % perspective of animal 2
    ma1 = allMOUSES{rec,2};
    ma2 = allMOUSES{rec,1};
    maS = cat(3,ma1,ma2);
    nn1 = size(maS,1);
    p1Dist = zeros(nx^2,size(maS,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(maS(:,:,Xi(i))),squeeze(maS(:,:,Yi(i))));
    end
    p1Dsmooth = zeros(size(p1Dist));
    for i = 1:size(p1Dist,1)
        p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
    end
    p1Dist = p1Dsmooth(IDX,:)';

    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs6;

    PJ = [d2anyt d1anyt];
    PJN = [d2anyN d1anyN];
    [W,~] = findWavelets(proj,9,parameters);
    data2 = log(W); data2(data2<-3) = -3;

    angInfo2 = zeros(size(angInfo));
    angInfo2(:,1:3) = angInfo(:,4:6); angInfo2(:,4:6) = angInfo(:,1:3);

    nnData = [data2 .25*rat2z .25*rat1z 2*p2v(:,[1 11 12 15 19 23]) 2*p1v(:,[1 11 12 15 19 23]) 10*angInfo2 .05*PJ 2*PJN];
    yData = tsne(nnData(1:20:end,:));
    amps = sum(nnData',1);
    dataMat2 = nnData;
    %         % FIND GOOD TEMPLATES AND SAVE THEM
    [signalData,signalAmps] = findTemplatesFromData(...
        nnData(1:20:end,:),yData,amps(1:20:end,:),numPerDataSet,parameters);
    mDC{rec,2} = signalData; mDCA{rec,2} = signalAmps;
end


% save('socSignalDataAmps_0520.mat','mDC','mDCA');

allD = combineCells(mDC(:),1); allA = combineCells(mDCA(:));
for i = 1:10
    tic
    Y{i} = tsne(allD);
    toc
end

for i = 1:10
    subplot(3,4,i);
    scatter(Y{i}(:,1),Y{i}(:,2),[],allD(:,end-12),'.');
    axis equal
end

yData = Y{9};
% save('train0520_soc.mat','yData','allD')

%% Reembed Soc
load('vecsValsSOCMouse0520.mat','C','L','mus','vals','vecs');
load('train0520_soc.mat','yData','allD')

vecs6 = vecs(:,1:6);
minF = .2; maxF = 5; pcaModes = 10; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = 400;
parameters.batchSize = 10000;

trainingSetData = allD; trainingEmbeddingZ = yData;

for rec = 1:length(allMOUSES)
    rec
    ma1 = allMOUSES{rec,1};
    ma2 = allMOUSES{rec,2};
    maS = cat(3,ma1,ma2);
    allz = squeeze(maS(:,3,[19 23 42 46])); zz = allz(:);
    fz = prctile(zz,10);

    nn1 = size(maS,1);
    p1Dist = zeros(nx^2,size(maS,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(maS(:,:,Xi(i))),squeeze(maS(:,:,Yi(i))));
    end
    p1Dsmooth = zeros(size(p1Dist));
    for i = 1:size(p1Dist,1)
        p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
    end
    p1Dist = p1Dsmooth(IDX,:)';

    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs6;

    % velocity information
    p1v = zeros(nn1,3); p2v = zeros(nn1,3);
    for i = 1:23
        p1v(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma1(:,:,i))).^2,2)),10)];
        p2v(:,i) = [0; medfilt1(sqrt(sum(diff(squeeze(ma2(:,:,i))).^2,2)),10)];
    end

    p1v(p1v>5) = 5; p2v(p2v>5) = 5;
    n1 = ma1(:,:,1); s1 = ma1(:,:,4); c1 = ma1(:,:,5); t1 = ma1(:,:,7);
    n2 = ma2(:,:,1); s2 = ma2(:,:,4); c2 = ma2(:,:,5); t2 = ma2(:,:,7);
    distc = sqrt(sum((c1-c2).^2,2));

    % angle information
    angInfo = zeros(nn1,6);
    vec1 = n2-s1; vec2 = n1-s1; angInfo(:,1) = findAngleDiff(vec1,vec2);
    vec1 = c2-s1; vec2 = n1-s1; angInfo(:,2) = findAngleDiff(vec1,vec2);
    vec1 = t2-s1; vec2 = n1-s1; angInfo(:,3) = findAngleDiff(vec1,vec2);
    vec1 = n1-s2; vec2 = n2-s2; angInfo(:,4) = findAngleDiff(vec1,vec2);
    vec1 = c1-s2; vec2 = n2-s2; angInfo(:,5) = findAngleDiff(vec1,vec2);
    vec1 = t1-s2; vec2 = n2-s2; angInfo(:,6) = findAngleDiff(vec1,vec2);
    angInfo = pi - angInfo;

    d2n = sqrt(sum((n1-n2).^2,2)); % dist of nose1 to other rat
    d2t = sqrt(sum((n1-t2).^2,2));
    d2c = sqrt(sum((n1-c2).^2,2));
    d2anyt = [d2n d2c d2t];
    d2anyN = d2anyt./distc;

    d1n = sqrt(sum((n2-n1).^2,2)); % dist of nose2 to other rat
    d1t = sqrt(sum((n2-t1).^2,2));
    d1c = sqrt(sum((n2-c1).^2,2));
    d1anyt = [d1n d1c d1t];
    d1anyN = d1anyt./distc;

    PJ = [d1anyt d2anyt];
    PJN = [d1anyN d2anyN];
    [W,~] = findWavelets(proj,9,parameters);
    data2 = log(W); data2(data2<-3) = -3;

    % compute other social features and make final feature matrix
    rat1z = squeeze(ma1(:,3,[1 4 5 8 11 12 15]))-fz;
    rat2z = squeeze(ma2(:,3,[1 4 5 8 11 12 15]))-fz;
    nnData = [data2 .25*rat1z .25*rat2z 2*p1v(:,[1 11 12 15 19 23]) 2*p2v(:,[1 11 12 15 19 23]) 10*angInfo .05*PJ 2*PJN];
    dataMat1 = nnData;

    % perspective of animal 2
    ma1 = allMOUSES{rec,2};
    ma2 = allMOUSES{rec,1};
    maS = cat(3,ma1,ma2);
    nn1 = size(maS,1);
    p1Dist = zeros(nx^2,size(maS,1));
    for i = 1:size(p1Dist,1)
        p1Dist(i,:) = returnDist3d(squeeze(maS(:,:,Xi(i))),squeeze(maS(:,:,Yi(i))));
    end
    p1Dsmooth = zeros(size(p1Dist));
    for i = 1:size(p1Dist,1)
        p1Dsmooth(i,:) = smooth(medfilt1(p1Dist(i,:),3),3);
    end
    p1Dist = p1Dsmooth(IDX,:)';

    p1 = bsxfun(@minus,p1Dist,mus);
    proj = p1*vecs6;

    PJ = [d2anyt d1anyt];
    PJN = [d2anyN d1anyN];
    [W,~] = findWavelets(proj,9,parameters);
    data2 = log(W); data2(data2<-3) = -3;

    angInfo2 = zeros(size(angInfo));
    angInfo2(:,1:3) = angInfo(:,4:6); angInfo2(:,4:6) = angInfo(:,1:3);

    nnData = [data2 .25*rat2z .25*rat1z 2*p2v(:,[1 11 12 15 19 23]) 2*p1v(:,[1 11 12 15 19 23]) 10*angInfo2 .05*PJ 2*PJN];
    dataMat2 = nnData;

    fprintf(1,'Finding Embeddings\n');
    [zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
        findTDistributedProjections_fmin(dataMat1,trainingSetData,...
        trainingEmbeddingZ,[],parameters);

    z1 = zValues; z1(~inConvHull,:) = zGuesses(~inConvHull,:);
    inCH1 = inConvHull;

    fprintf(1,'Finding Embeddings\n');
    [zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
        findTDistributedProjections_fmin(dataMat2,trainingSetData,...
        trainingEmbeddingZ,[],parameters);

    z2 = zValues; z2(~inConvHull,:) = zGuesses(~inConvHull,:);
    inCH2 = inConvHull;

    filename = allMOUSES_files{rec}; 
    % save(['RE_soc/RE_SOC_' num2str(rec) '.mat'],'z1','inCH1','z2','inCH2','filename')
end

%%
% make lone LL, WR, diff
% allMOUSEL = [Bpred; Wpred];
% allMOUSES = [BBpred; WWpred; BWpred];
load('/Users/ugne/Library/CloudStorage/Dropbox/MultiFlyAnalysis_OLD/visualization/colormaps.mat')

EVL = cell(128,1); ICVL = zeros(128,1);
for i = 1:128
    i
    load(['RE_lone/RE_LONE_' num2str(i) '.mat'],'z','inConvHull','filename')
    EVL{i} = z;
    ICVL(i) = mean(inConvHull);
end

evall = combineCells(EVL);
[xx d] = findPointDensity(evall,1,501,[-65 65]);
D = d;
 % Watershed
LL = watershed(-d,8);
LL2 = LL; LL2(d < 1e-6) = -1;
LL3 = zeros(size(LL));
for i = 1:501
    for j = 1:501
        if LL2(i,j)==0
            LL3(i,j)=1;
        end
    end
end
LLBW = LL2==0;
LLBWB = bwboundaries(LLBW);
llbwb = LLBWB(2:end);
llbwb = combineCells(llbwb');
figure(4); imagesc(D); axis equal off; colormap(flipud(gray)); caxis([0 6e-4]);
hold on; scatter(llbwb(:,2),llbwb(:,1),'.','k'); % loneMapGray

vSmooth = .5;
medianLength = 1;
pThreshold = [];
minRest = 5;
obj = [];
fitOnly = true;
numGMM = 2;

[wr,segments,v,obj,pRest,vals,vx,vy] = ...
    findWatershedRegions_v2(evall,xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);

load('reembeddingInfoLone.mat','xx','LL','vSmooth','medianLength','pThreshold','minRest','obj','numGMM')

WRL = cell(size(EVL));
for i = 1:128
    [WRL{i},segments,v,obj,pRest,vals,vx,vy] = ...
        findWatershedRegions_v2(EVL{i},xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,false,numGMM);
end


sc = scM;
offx = 150;
offy = 0;
offz = 150;
coffx = zeros(4,5); coffz = zeros(4,5);
for i = 1:5
    for j = 1:4
        coffx(i,j) = (i-1)*offx;
        coffz(i,j) = (j-1)*offz;
    end
end

[groups,~,~] = makeGroupsAndSegments(WRL(:),max(max(LL)),ones(1,128),15);
allMOUSE_A = cell(size(allMOUSE));
for i = 1:length(allMOUSE)
    [allMOUSE_A{i},~] = alignDannceNF(allMOUSE{i});
end

M = allMOUSE_A;

for tm = 1:max(max(LL))
    try
        G = groups{tm};
        lgs = size(G,1);
        if lgs>20
            newG = G(randperm(lgs,20),:);
        else
            newG = G;
        end
        lenG = 60;
        gMarkers = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = M{newG(i,1)}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers{i} = tgm(1:lenG,:,:);
            end
        end
        allM = []; allJ = []; allC = []; allmc = [];
        for i = 1:size(newG,1)
            cx = coffx(i); cz = coffz(i);
            NM = gMarkers{i};
            NM(:,1,:) = NM(:,1,:)+cx;
            NM(:,3,:) = NM(:,3,:)+cz;
            
            allM = cat(3,allM,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ = cat(1,allJ,newJ);
            allmc = cat(1,allmc,sc);
            allC = cat(1,allC,zeros(23,3));
        end
        sk2.color = allC; sk2.joints_idx = allJ; sk2.mcolor = allmc;
        
        close all;
        findic = figure('Name','Rat Test');
        h = cell(1,1);
        h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'MarkerSize',30,'LineWidth',2);
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition', [0.5174 0.1267 0.2398 0.7226]);
        view(h{1}.getAxes,20,30); axis equal
        axis([-50 650 -50 50 -100 550]);
        set(gca,'Color','w'); %set(gca,'Color','k')
        set(gcf,'Color','white'); %set(gcf,'Color','black');
        
        savePath = ['bradies_lone/lone_a_' num2str(tm) '.avi'];
        frames = 2:lenG;
        
        % Uncomment to write the Animation to video.
        % h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
    catch
    end
end

M = allMOUSE;
for tm = 1:max(max(LL))
    try
        G = groups{tm};
        lgs = size(G,1);
        if lgs>20
            newG = G(randperm(lgs,20),:);
        else
            newG = G;
        end
        lenG = 60;
        gMarkers = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = M{newG(i,1)}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers{i} = tgm(1:lenG,:,:);
            end
        end
        allM = []; allJ = []; allC = []; allmc = [];
        for i = 1:size(newG,1)
            cx = 0; cz = 0; %cx = coffx(i); cz = coffz(i);
            NM = gMarkers{i};
            NM(:,1,:) = NM(:,1,:)+cx;
            NM(:,3,:) = NM(:,3,:)+cz;
            
            allM = cat(3,allM,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ = cat(1,allJ,newJ);
            allmc = cat(1,allmc,sc);
            allC = cat(1,allC,zeros(23,3));
        end
        sk2.color = allC; sk2.joints_idx = allJ; sk2.mcolor = allmc;
        
        close all;
        findic = figure('Name','Rat Test');
        h = cell(1,1);
        h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'MarkerSize',30,'LineWidth',2);
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition', [0.5174 0.1267 0.2398 0.7226]);
        view(h{1}.getAxes,20,30); axis equal
        axis([-160 250 -100 300 0 125]);
        set(gca,'Color','w'); %set(gca,'Color','k')
        set(gcf,'Color','white'); %set(gcf,'Color','black');
        
        savePath = ['bradies_lone/lone_s_' num2str(tm) '.avi'];
        frames = 2:lenG;
        
        % Uncomment to write the Animation to video.
        % h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
    catch
    end
end





%%
% make LLsoc, WRsoc, diff

EVS = cell(48,2); ICVS = zeros(48,2);
for i = 1:48
    i
    load(['RE_soc/RE_SOC_' num2str(i) '.mat'],'z1','z2','inCH1','inCH2','filename');
    EVS{i,1} = z1; EVS{i,2} = z2;
    ICVS(i,1) = mean(inCH1); ICVS(i,2) = mean(inCH2);
end

evallsoc = combineCells(EVS(:));
[xxsoc, dsoc] = findPointDensity(evallsoc,1.5,501,[-85 85]);
 % Watershed
LLsoc = watershed(-dsoc,8);
LL2soc = LLsoc; LL2soc(dsoc < 1e-6) = -1;
LL3soc = zeros(size(LLsoc));
for i = 1:501
    for j = 1:501
        if LL2soc(i,j)==0
            LL3soc(i,j)=1;
        end
    end
end
LLBWsoc = LL2soc==0;
LLBWBsoc = bwboundaries(LLBWsoc);
llbwbsoc = LLBWBsoc(2:end);
llbwbsoc = combineCells(llbwbsoc');
figure(4); imagesc(dsoc); axis equal off; colormap(cmap1); caxis([0 3e-4]);
hold on; scatter(llbwbsoc(:,2),llbwbsoc(:,1),1,'k');

vSmooth = .5;
medianLength = 1;
pThreshold = [];
minRest = 10;
obj = [];
fitOnly = true;
numGMM = 2;

[wrsoc,segments,v,obj,pRest,vals,vx,vy] = ...
    findWatershedRegions_v2(evallsoc,xxsoc,LL2soc,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);

WRS = cell(size(EVS));
for i = 1:48
    for j = 1:2
        [WRS{i,j},segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(EVS{i,j},xxsoc,LL2soc,vSmooth,medianLength,pThreshold,minRest,obj,false,numGMM);
    end
end


[groupsS,~,~] = makeGroupsAndSegments(WRS(:),max(max(LL2soc)),ones(96,1),15);
sc = scM;
offx = 450;
offy = 450;
offz = 0;
coffx = zeros(3,3); coffy = zeros(3,3); coffz = zeros(3,3);
for i = 1:3
    for j = 1:3
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy';



p2 = [allMOUSES(:,2) allMOUSES(:,1)];
PRED = [allMOUSES; p2];


tm = 1;
isrr = ones(1,96);
for tm = 1:133
    try
        G = groupsS{tm};
        lgs = size(G,1);
        if lgs>9
            newG = G(randperm(lgs,9),:);
        else
            newG = G;
        end
        lenG = 150;
        gMarkers1 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = PRED{newG(i,1),1}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers1{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers1{i} = tgm(1:lenG,:,:);
            end
        end
        gMarkers2 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = PRED{newG(i,1),2}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers2{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers2{i} = tgm(1:lenG,:,:);
            end
        end

        allM = []; allJ = []; allC = []; allmc = [];
        for i = 1:size(newG,1)
            cx = coffx(i); cy = coffy(i); cz = coffz(i);
            NM = gMarkers1{i};
            NM(:,1,:) = NM(:,1,:)+cx;
            NM(:,2,:) = NM(:,2,:)+cy;
            NM(:,3,:) = NM(:,3,:)+cz;
            allM = cat(3,allM,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ = cat(1,allJ,newJ);
            allC = cat(1,allC,skeleton.color);
        end
        sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));

        allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
        for i = 1:size(newG,1)
            cx = coffx(i); cy = coffy(i); cz = coffz(i);
            NM = gMarkers2{i};
            NM(:,1,:) = NM(:,1,:)+cx;
            NM(:,2,:) = NM(:,2,:)+cy;
            NM(:,3,:) = NM(:,3,:)+cz;
            allM2 = cat(3,allM2,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ2 = cat(1,allJ,newJ);
            allC2 = cat(1,allC2,skeleton.color);
        end
        sk3.mcolor = zeros(size(allC2)); sk3.joints_idx = allJ; sk3.color = zeros(size(allC));

        close all;
        findic = figure('Name','Rat Test');
        h = cell(2,1);
        h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
        h{2} = Keypoint3DAnimator(allM2,sk3,'Axes',h{1}.Axes,'lineWidth',1.5);
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition',[0.3934 0.0679 0.5402 0.8232]);
                view(h{1}.getAxes,20,30); axis equal

        x1 = 45; y1 = 90;
        for i = 1:3
            x = x1+coffx(i,1);
            for j = 1:3
                y = y1+coffy(1,j);
                r  = 215;
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + x;
                yunit = r * sin(th) + y;
                h2 = plot(xunit, yunit,'Color',[.1 .5 .1]); hold on;
            end
        end
        axis([-200 1205 -200 1205 -10 100]);
        set(gcf,'Color','white');

        savePath = ['bradies_soc2/soc_s_' num2str(tm) '.avi'];
        frames = 1:lenG;

        % Uncomment to write the Animation to video.
        % h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
        pause(2)
    catch
    end
end



%%
wrFINE = cell(size(WRL));
for i = 1:length(wrFINE)
    wrt = WRL{i};
    wrt(wrt==0) = nan; 
    wrt = fillmissing(wrt,'nearest');
    wrFINE{i} = wrt;
end

wrsocFINE = cell(size(WRS));
for i = 1:length(wrsocFINE)
    for j = 1:2
        wrt = WRS{i,j}; wrt(wrt==134) = 1;
        wrt(wrt==0) = nan;
        wrt = fillmissing(wrt,'nearest');
        wrsocFINE{i,j} = wrt;
    end
end

wrCOARSE = cell(size(wrFINE));
for i = 1:length(wrCOARSE)
   wrt = wrFINE{i};
    wrCOARSE{i} = labelsL(wrt);
end

wrsocCOARSE = cell(size(wrsocFINE));
for i = 1:length(wrsocCOARSE)
    for j = 1:2
        wrt = wrsocFINE{i,j}; wrt(wrt==134) = 1;
        wrsocCOARSE{i,j} = labelsS(wrt);
    end
end

% B/W, mouseID, partner B/W, partnermouseID
mouseIDL = cell(32,4);
mouseIDS = cell(48,4);
for i = 1:32
    t = allMOUSEL_files{i};
    mouseIDL{i,1} = t(end-1);
    mouseIDL{i,2} = t(end);
end

for i = 1:48
    t = allMOUSES_files{i};
    mouseIDS{i,1} = t(end-4);
    mouseIDS{i,2} = t(end-3);
    mouseIDS{i,3} = t(end-1);
    mouseIDS{i,4} = t(end);
end


for i = 1:48
    f = allMOUSES_files{i};
    t1 = f(end-4:end-3); t2 = f(end-1:end);
    mouseTagsS{i,1} = t1; mouseTagsS{i,2} = t2;
end
load('mouse_embedding_info.mat','wrFINE','wrsocFINE','LL','LL2','LLC','llbwb','LLsoc','LL2soc','LLCsoc','llbwbsoc',...
    'mouseOrderL','mouseOrderS','allMOUSEL_files','allMOUSES_files','mouseIDL','mouseIDS');

load('mouse_embedding_info_2.mat','loneUsage','mouseIDLone','mouseIDSoc','mouseOrderL','mouseIDall','loneUsageS',...
    'mIDSoc','socUsage1','socUsage2')



llt = 30000;
maxL = max(max(LL)); maxS = max(max(LL2soc));
loneUsage = zeros(128,maxL); socUsage1 = zeros(48,maxS); socUsage2 = zeros(48,maxS);

for i = 1:128
    tt = hist(wrFINE{i},1:double(maxL))./llt;
    loneUsage(i,:) = tt;
end

for i = 1:48
    tt = hist(wrsocFINE{i,1},1:double(maxS))./llt;
    socUsage1(i,:) = tt;
    tt2 = hist(wrsocFINE{i,2},1:double(maxS))./llt;
    socUsage2(i,:) = tt2;
end

socL1 = wrFINE(33:80); socL2 = wrFINE(81:128);
SL1 = zeros(48,maxL); SL2 = zeros(48,maxL);
for i = 1:48
    tt = hist(socL1{i},1:double(maxL))./llt;
    SL1(i,:) = tt;
    tt2 = hist(socL2{i},1:double(maxL))./llt;
    SL2(i,:) = tt2;
end


loneUsageC = zeros(128,14); socUsage1C = zeros(48,7); socUsage2C = zeros(48,7);

for i = 1:128
    tt = hist(wrCOARSE{i},1:double(14))./llt;
    loneUsageC(i,:) = tt;
end

for i = 1:48
    tt = hist(wrsocCOARSE{i,1},1:7)./llt;
    socUsage1C(i,:) = tt;
    tt2 = hist(wrsocCOARSE{i,2},1:7)./llt;
    socUsage2C(i,:) = tt2;
end

socL1C = wrCOARSE(33:80); socL2C = wrCOARSE(81:128);
SL1C = zeros(48,14); SL2C = zeros(48,14);
for i = 1:48
    tt = hist(socL1C{i},1:14)./llt;
    SL1C(i,:) = tt;
    tt2 = hist(socL2C{i},1:14)./llt;
    SL2C(i,:) = tt2;
end








EVS1 = zeros(48,30000); EVS2 = zeros(48,30000);
for i = 1:48
    t = socL1{i}; tt = labelsL(t);
    EVS1(i,:) = tt;
    t2 = socL2{i}; tt2 = labelsL(t2);
    EVS2(i,:) = tt2;
end
figure(1); imagesc(EVS1); hold on; colormap(colorsCoarse); caxis([0 14]) 
figure(2); imagesc(EVS2); hold on; colormap(colorsCoarse); caxis([0 14]) 






EVL = zeros(32,30000);
for i = 1:32
    t = wrFINE{i}; tt = labelsL(t);
    EVL(i,:) = tt;
end
imagesc(EVL); hold on; colormap(colorsCoarse); caxis([0 14]) 



mouseDC = cell(48,1); 
for i = 1:48
    ma1 = allMOUSES{i,1}; ma2 = allMOUSES{i,2};
    p1_xy = ma1(:,:,[1 5 7]);
    p2_xy = ma2(:,:,[1 5 7]);
    n1 = p1_xy(:,:,1); c1 = p1_xy(:,:,2); t1 = p1_xy(:,:,3);
    n2 = p2_xy(:,:,1); c2 = p2_xy(:,:,2); t2 = p2_xy(:,:,3);
    distc = sqrt(sum((c1-c2).^2,2));
    mouseDC{i} = distc;
end

mouseCV = cell(128,1);
for i = 1:128
    pred = allMOUSE{i};
    jd = [0; smooth(sqrt(sum(diff(squeeze(pred(:,:,7))).^2,2)),20)];
    mouseCV{i} = jd;
end

hcv = zeros(128,20);
for i = 1:128
    hh = hist(mouseCV{i},.25:.5:10);
    hcv(i,:) = hh./(sum(hh));
end


mouseTagsS = cell(48,2);
mouseTagsL = cell(32,1);

for i = 1:32
    f = allMOUSEL_files{i};
    t1 = f(end-1:end);
    mouseTagsL{i} = t1;
end
for i = 1:48
    f = allMOUSES_files{i};
    t1 = f(end-4:end-3); t2 = f(end-1:end);
    mouseTagsS{i,1} = t1; mouseTagsS{i,2} = t2;
end

isW_L = zeros(size(mouseTagsL));
isW_S = zeros(size(mouseTagsS));

for i = 1:32
    if strfind(mouseTagsL{i},'W')
        isW_L(i) = 1;
    end
end
for i = 1:48
    if strfind(mouseTagsS{i,1},'W')
        isW_S(i,1) = 1;
    end
    if strfind(mouseTagsS{i,2},'W')
        isW_S(i,2) = 1;
    end
end

bb_rec = find(isW_S(:,1)==0 & isW_S(:,2)==0);
ww_rec = find(isW_S(:,1)==1 & isW_S(:,2)==1);
bw_rec = find(isW_S(:,1)==0 & isW_S(:,2)==1);

bbdc = combineCells(mouseDC(bb_rec));
wwdc = combineCells(mouseDC(ww_rec));
bwdc = combineCells(mouseDC(bw_rec));

hbb = hist(bbdc,0:10:400);
hww = hist(wwdc,0:10:400);
hbw = hist(bwdc,0:10:400);

allDC = [combineCells(mouseDC); combineCells(mouseDC)];
allSB = [combineCells(wrsocFINE(:,1)); combineCells(wrsocFINE(:,2))];

socDist = zeros(1,133);
for i = 1:133
    idx = find(allSB==i);
    mdc = median(allDC(idx),'omitnan');
    socDist(i) = mdc;
end

plot(hbb./sum(hbb),'LineWidth',3,'Color','k'); hold on;
plot(hww./sum(hww),'LineWidth',3,'Color','r'); hold on;
plot(hbw./sum(hbw),'LineWidth',3,'Color',[.5 .5 .5]); 

%%
idle = [29 31 36 39];
slow = [44 58 62 66 68 78 85 92];
head = [74 79 82];
groom = [33 37 40 46 52];
crouched = [15 16 21 23 24 28];
actCrouch = [20 22 25 26 27 30 32 34 35 38 41 43 47 48 50 56];
stepsCrouched = [42 45 49 51 53 54 55 59 61 ];
rear = [13 14 17 18 19];
highRear = [1 2 3 4 5 6 7 8 9 10 11 12]; 
slowExp = [63 76 88 90 91];
Exp = [57 60 64 72 73 75 80 81 87 93 97 99 100 103 104 106 107 111 114];
stepExp = [69 83 86 101 108 110 112 113 115 117 118]; 
locSlow = [67 71 77 84 94 95 96 98 105];
locFast = [65 70 89 102 109 116];

labelsL = zeros(1,118);
labelsL(idle) = 1;
labelsL(slow) = 2;
labelsL(head) = 3;
labelsL(groom) = 4;
labelsL(crouched) = 5;
labelsL(actCrouch) = 6;
labelsL(stepsCrouched) = 7;
labelsL(rear) = 8;
labelsL(highRear) = 9;
labelsL(slowExp) = 10;
labelsL(Exp) = 11;
labelsL(stepExp) = 12;
labelsL(locSlow) = 13;
labelsL(locFast) = 14;

%% diff Black-White, Lone-social B/W, 

    colorsCoarse = [1 1 1;
        0 .3333 1;
        0 .7222 1;
        .1333 1 .8667;
        .6 1 .4;
        1 .9375 0;
        1 .6 0;
        1 0 0;
        .6 0 0;
        .4 0 0];
    

    colorsCoarse = [1 1 1;
        0 .16 1; % idle
        0 .3833 1; %slow
        0 .7222 1; % head
        .1333 1 .8667; %groom
        .6 1 .4; % crouched
        1 .9375 0; % act crouchd
        1 .75 0; % steps crouched
        1 .6 0; % reared
        1 .4 0; % high rear
        1 0 .2; % slow exp
        1 0 .4; % explore
        1 0 0; % step explore
        .6 0 0; % loc Slow
        .4 0 0]; % loc Fast
    

LLC = zeros(size(LL2));
for i = 1:118
    cid = find(LL2==i);
    LLC(cid) = labelsL(i);
end

imagesc(LLC); hold on; colormap(cmapLL); caxis([0 15])
for i = 1:118
    [xi yi] = find(LL2==i);
    mxi = mean(xi); myi = mean(yi);
    text(myi,mxi,num2str(i));
end

allcv = combineCells(mouseCV);
allwr = combineCells(wrFINE);
meancv = zeros(118,1);

for b = 1:118
    idx = find(allwr==b);
    mdc = median(allcv(idx),'omitnan');
    meancv(b) = mdc;
end

mouseOrderL = []; mouseOrderLTag = [];
for b = 1:14
    xi = find(labelsL==b);
    dxi = meancv(xi);
    [dsorted, dsxi] = sort(dxi,'ascend');
    sortedxi = xi(dsxi)';
    mouseOrderL = [mouseOrderL; sortedxi];
    mouseOrderLTag = [mouseOrderLTag; b*ones(size(sortedxi))];
end


b = bar(meancv(mouseOrderL));
b.FaceColor = 'flat';
for i = 1:length(mouseOrderL)
    b.CData(i,:) = colorsCoarse(mouseOrderLTag(i)+1,:);
end

imagesc(loneUsage(:,mouseOrderL))


mouseBLd1 = loneUsage(1:8,mouseOrderL);
mouseBLd2 = loneUsage(9:16,mouseOrderL);
mouseWLd1 = loneUsage(17:24,mouseOrderL);
mouseWLd2 = loneUsage(25:32,mouseOrderL);

imagesc(log([mouseBLd1; inf(1,118); mouseBLd2; inf(3,118); ...
    mouseWLd1; inf(1,118); mouseWLd2]));
colormap(gray);


SL1_spaced = [SL1(1:12,mouseOrderL); inf(1,118); SL1(13:28,mouseOrderL); inf(1,118); SL1(29:48,mouseOrderL)];
imagesc(log(SL1_spaced)); colormap(gray)

SL2_spaced = [SL2(1:12,mouseOrderL); inf(1,118); SL2(13:28,mouseOrderL); inf(1,118); SL2(29:48,mouseOrderL)];
imagesc(log(SL2_spaced)); colormap(gray)

SS1_spaced = [socUsage1(1:12,mouseOrderS); inf(1,133); socUsage1(13:28,mouseOrderS); inf(1,133); socUsage1(29:48,mouseOrderS)];
imagesc(log(SS1_spaced)); colormap(gray)

SS2_spaced = [socUsage2(1:12,mouseOrderS); inf(1,133); socUsage2(13:28,mouseOrderS); inf(1,133); socUsage2(29:48,mouseOrderS)];
imagesc(log(SS2_spaced)); colormap(gray)



% idle, slow, head,  groom, crouched, actCrouched, stepsCrouched, rear,
% highRear, slowExplore, explore, step explore, locSlow, locFast



imagesc(LLC); hold on; colormap(colorsCoarse); caxis([0 14])
hold on; scatter(llbwb(:,2),llbwb(:,1),'.','k'); axis equal; % LLClone.tif

imagesc(D); colormap(flipud(gray)); caxis([0 6e-4]);
hold on; scatter(llbwb(:,2),llbwb(:,1),'.','k'); axis equal % loneMapGray


sNI = [41 42 43 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 ...
    62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 ...
    81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 ...
    97 98 99 101 102 103 104 105 108 109 110 111 ...
    113 114 115 116 118 119 120 122 124 125 127 128 129 130 132];

sLM1 = [20 24 28 30 32 36 37 38 39 40]; 
sLM2 = [31 89 100 106 107 133];

sM1 = [17 19 22 25 26 27 29 33 34 44];
sM2 = [112 117 121 123 126 131];
sLMI = [31 35 43];
sMI = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 21 23];

mutRear = [12 15 19 92 94 101 103 105 116];

labelsS = zeros(1,133);
labelsS(sNI) = 1;
labelsS(sLM1) = 2;
labelsS(sM1) = 3;
labelsS(sLM2) = 4;
labelsS(sM2) = 5;
labelsS(sLMI) = 6;
labelsS(sMI) = 7;


LLCsoc = zeros(size(LL2soc));
for i = 1:133
    cid = find(LL2soc==i);
    LLCsoc(cid) = labelsS(i);
end

cmapSOC = [1 1 1;
    .8 .8 .8;
     0 0 .7;
     0 0 1;
     .7 0 0; 
     1 0 0;
     .6 0 .6;
     .8 0 .8];

figure(4); imagesc(dsoc); axis equal off; colormap(cmap1); caxis([0 3e-4]);
hold on; scatter(llbwbsoc(:,2),llbwbsoc(:,1),1,'k','.');
set(gcf,'Position',[653 242 868 759]);

imagesc(LLCsoc); hold on; colormap(cmapSOC); caxis([0 8])
hold on; scatter(llbwbsoc(:,2),llbwbsoc(:,1),1,'k','.'); axis equal off
set(gcf,'Position',[653 242 868 759]);


for i = 1:133
    [xi yi] = find(LL2soc==i);
    mxi = mean(xi); myi = mean(yi);
    text(myi,mxi,num2str(i));
end

%socDist


mouseOrderS = []; mouseOrderSTag = [];
for b = 1:7
    xi = find(labelsS==b);
    dxi = socDist(xi);
    [dsorted, dsxi] = sort(dxi,'descend');
    sortedxi = xi(dsxi)';
    mouseOrderS = [mouseOrderS; sortedxi];
    mouseOrderSTag = [mouseOrderSTag; b*ones(size(sortedxi))];
end
mouseOrderS2 = mouseOrderS;

b = bar(socDist(mouseOrderS));
b.FaceColor = 'flat';
for i = 1:length(mouseOrderS)
    b.CData(i,:) = cmapSOC(mouseOrderSTag(i)+1,:);
end

%% mouse stats

load('mouse_embedding_info.mat','wrFINE','wrsocFINE','LL','LL2','LLC','llbwb','LLsoc','LL2soc','LLCsoc','llbwbsoc',...
    'mouseOrderL','mouseOrderS','allMOUSEL_files','allMOUSES_files','mouseIDL','mouseIDS');

mouseIDLone = cell(size(allMOUSEL_files));
for i = 1:length(mouseIDLone)
    mouseIDLone{i} = allMOUSEL_files{i}(end-1:end);
end
mouseIDSoc = cell(size(allMOUSES_files,1),2);
for i = 1:length(mouseIDSoc)
    mouseIDSoc{i,1} = allMOUSES_files{i}(end-4:end-3);
    mouseIDSoc{i,2} = allMOUSES_files{i}(end-1:end);
end

mouseBLd1 = loneUsage(1:8,mouseOrderL);
mouseBLd2 = loneUsage(9:16,mouseOrderL);
mouseWLd1 = loneUsage(17:24,mouseOrderL);
mouseWLd2 = loneUsage(25:32,mouseOrderL);


Blone = loneUsage(1:16,mouseOrderL);
BloneID = mouseIDLone(1:16);
Wlone = loneUsage(17:32,mouseOrderL);
WloneID = mouseIDLone(17:32);

% BB - 1:12, WW - 16?
loneUsageS = loneUsage(33:end,mouseOrderL); loneUsageSC = loneUsageC(33:end,:);
mouseIDall = [mouseIDLone; mouseIDSoc(:)];
mIDSoc = mouseIDall(33:end);

BBsl = loneUsageS([1:12 49:60],:);
WWsl = loneUsageS([13:28 61:76],:);
BWsl = loneUsageS(29:48,:);
WBsl = loneUsageS(77:96,:);


cBBsl = loneUsageSC([1:12 49:60],:);
cWWsl = loneUsageSC([13:28 61:76],:);
cBWsl = loneUsageSC(29:48,:);
cWBsl = loneUsageSC(77:96,:);


BBslid = mIDSoc([1:12 49:60]);
WWslid = mIDSoc([13:28 61:76]);
BWslid = mIDSoc(29:48);
WBslid = mIDSoc(77:96);


[pvalue, testStat, lme] = testGroupLme(Blone, Wlone, BloneID, WloneID); 
imagesc(testStat); colormap(cmap1); caxis([-.01 .22]); hold on; axis off equal


LMpvals = nan(4,4); LMstats = nan(4,4);
[pvalue, testStat, lme] = testGroupLme(BBsl, WWsl, BBslid, WWslid); LMpvals(1,2) = pvalue; LMstats(1,2) = testStat; %6.1281e-13
[pvalue, testStat, lme] = testGroupLme(BBsl, BWsl, BBslid, BWslid); LMpvals(1,3) = pvalue; LMstats(1,3) = testStat; %0.0115
[pvalue, testStat, lme] = testGroupLme(BBsl, WBsl, BBslid, WBslid); LMpvals(1,4) = pvalue; LMstats(1,4) = testStat; %7.9680e-15

[pvalue, testStat, lme] = testGroupLme(WWsl, BWsl, WWslid, BWslid); LMpvals(2,3) = pvalue; LMstats(2,3) = testStat; %3.0941e-16
[pvalue, testStat, lme] = testGroupLme(WWsl, WBsl, WWslid, WBslid); LMpvals(2,4) = pvalue; LMstats(2,4) = testStat; %0.2994

[pvalue, testStat, lme] = testGroupLme(WBsl, BWsl, WBslid, BWslid); LMpvals(3,4) = pvalue; LMstats(3,4) = testStat; %5.6237e-11

imagesc(LMstats); colormap(cmap1); caxis([-.01 .22]); hold on; axis off equal

for i = 1:4
    for j = 1:4
        if ~isnan(LMstats(i,j))
        text(j,i,num2str(LMstats(i,j)));
        end
    end
end



socUsageAll = [socUsage1(:,mouseOrderS); socUsage2(:,mouseOrderS)];
socUsageAllC = [socUsage1C; socUsage2C];


BBs = socUsageAll([1:12 49:60],:);
WWs = socUsageAll([13:28 61:76],:);
BWs = socUsageAll(29:48,:);
WBs = socUsageAll(77:96,:);

cBBs = socUsageAllC([1:12 49:60],:);
cWWs = socUsageAllC([13:28 61:76],:);
cBWs = socUsageAllC(29:48,:);
cWBs = socUsageAllC(77:96,:);

SMpvals = nan(4,4); SMstats = nan(4,4);

[pvalue, testStat, lme] = testGroupLme(BBs, WWs, BBslid, WWslid); SMpvals(1,2) = pvalue; SMstats(1,2) = testStat; %4.9600e-13
[pvals, fdrm] = findSigBeh(BBs,WWs,BBslid,WWslid);
diffB = mean(WWs)-mean(BBs);
up = mouseOrderSTag(find(fdrm<.05 & diffB>50/30000));
down = mouseOrderSTag(find(fdrm<.05 & diffB<-50/30000));
hup = hist(up,1:7); hdown = hist(down,1:);



[pvalue, testStat, lme] = testGroupLme(BBs, BWs, BBslid, BWslid); SMpvals(1,3) = pvalue; SMstats(1,3) = testStat; %1.5211e-14
[pvalue, testStat, lme] = testGroupLme(BBs, WBs, BBslid, WBslid); SMpvals(1,4) = pvalue; SMstats(1,4) = testStat; %2.7955e-14

[pvalue, testStat, lme] = testGroupLme(WWs, BWs, WWslid, BWslid); SMpvals(2,3) = pvalue; SMstats(2,3) = testStat; %8.7059e-09
[pvalue, testStat, lme] = testGroupLme(WWs, WBs, WWslid, WBslid); SMpvals(2,4) = pvalue; SMstats(2,4) = testStat; %5.1174e-09

[pvalue, testStat, lme] = testGroupLme(WBs, BWs, WBslid, BWslid); SMpvals(3,4) = pvalue; SMstats(3,4) = testStat; %2.0781e-12


% BB compared with WW
% BW compared with BB
% WB compared with WW

diffBB_WWsl = mean(WWsl)-mean(BBsl);
diffBW_BBsl = mean(BWsl)-mean(BBsl);
diffWB_WWsl = mean(WBsl)-mean(WWsl);
diffWB_BBsl = mean(WBsl)-mean(BBsl);



diffBB_WW = mean(WWs)-mean(BBs);
diffBW_BB = mean(BWs)-mean(BBs);
diffWB_WW = mean(WBs)-mean(WWs);

[x,wwbbsig] = findSigBeh(WWs,BBs,WWslid,BBslid);
[x,bwbbsig] = findSigBeh(BWs,BBs,BWslid,BBslid);
[x,wbwwsig] = findSigBeh(WBs,WWs,WBslid,WWslid);
ssig = [wwbbsig; bwbbsig; wbwwsig];

[x,wwbbsigl] = findSigBeh(WWsl,BBsl,WWslid,BBslid);
[x,bwbbsigl] = findSigBeh(BWsl,BBsl,BWslid,BBslid);
[x,wbwwsigl] = findSigBeh(WBsl,WWsl,WBslid,WWslid);
slsig = [wwbbsigl; bwbbsigl; wbwwsigl];

figure(1);
diffss = [diffBB_WW; diffBW_BB; diffWB_WW];
imagesc(diffss); colormap(cmapdiff); caxis([-.05 .05]); hold on;
for g = 1:3
    for b = 1:133
        if abs(diffss(g,b)) > 100/30000 && ssig(g,b)<.05
            scatter(b,g,30,'w','filled');
            scatter(b,g,15,'k','filled');
        end
    end
end
set(gcf,'Position',[45 835 1332 55]); axis off


figure(2); 
diffsl = [diffBB_WWsl; diffBW_BBsl; diffWB_WWsl];
imagesc(diffsl); colormap(cmapdiff); caxis([-.05 .05]); hold on;
for g = 1:3
    for b = 1:118
        if abs(diffsl(g,b)) > 100/30000 && slsig(g,b)<.05
            scatter(b,g,30,'w','filled');
            scatter(b,g,15,'k','filled');
        end
    end
end
set(gcf,'Position',[45 835 1332 55]); axis off

%% individual behaviors - mouse plots
% 105, 34, 43, 9

% 
% 
b = 23; %105 - far apart, mutual rear (black/black)
C = cell(4,1);
C{1} = BBs(:,b);
C{2} = BWs(:,b);
C{3} = WWs(:,b);
C{4} = WBs(:,b);
close all;
cents = [1 3 5 7];
ccol = [0 0 0; .7 .7 .7; 1 0 0; 0.8555 0.4883 0.8242];
[outinfo,tt] = mySwarmNA(C,cents,1,ccol);
axis([0 8 0 .07])
set(gca, 'Ticklength', [0 0])
set(gcf,'Position',[640 670 138 331]);


b = 92; %34 - aproach from behind (WB) - pink-dom
close all;
C = cell(4,1);
C{1} = BBs(:,b);
C{2} = BWs(:,b);
C{3} = WWs(:,b);
C{4} = WBs(:,b);
cents = [1 3 5 7];
[outinfo,tt] = mySwarmNA(C,cents,1,ccol);
axis([0 8 0 .08])
set(gca, 'Ticklength', [0 0])
set(gcf,'Position',[640 670 138 331]);


b = 112; %43
close all;
C = cell(4,1);
C{1} = BBs(:,b);
C{2} = BWs(:,b);
C{3} = WWs(:,b);
C{4} = WBs(:,b);
cents = [1 3 5 7];
[outinfo,tt] = mySwarmNA(C,cents,1,ccol);
axis([0 8 0 .06])
set(gca, 'Ticklength', [0 0])
set(gcf,'Position',[640 670 138 331]);


b = 129; %9 - mutual close - red dom (WW)
close all;
C = cell(4,1);
C{1} = BBs(:,b);
C{2} = BWs(:,b);
C{3} = WWs(:,b);
C{4} = WBs(:,b);
cents = [1 3 5 7];
[outinfo,tt] = mySwarmNA(C,cents,1,ccol);
axis([0 8 0 .07])
set(gca, 'Ticklength', [0 0])
set(gcf,'Position',[640 670 138 331]);




BBslid = mIDSoc([1:12 49:60]);
WWslid = mIDSoc([13:28 61:76]);
BWslid = mIDSoc(29:48);
WBslid = mIDSoc(77:96);

imagesc(SMstats); colormap(cmap1); caxis([-.01 .22]); hold on; axis off equal

for i = 1:4
    for j = 1:4
        if ~isnan(SMstats(i,j))
        text(j,i,num2str(LMstats(i,j)));
        end
    end
end


% 3.84 million mouse postures


%%  mouse social
% PCA for lone
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca([Blone; Wlone]); % 47.2224, 24.2416, 9.5787
colorLone = [repmat([0 0 0],[16,1]); repmat([1 0 0],[16,1])];
scatter(SCORE(:,1),SCORE(:,2),100,colorLone,'filled')
axis equal; 


BBsl = loneUsageS([1:12 49:60],:);
WWsl = loneUsageS([13:28 61:76],:);
BWsl = loneUsageS(29:48,:);
WBsl = loneUsageS(77:96,:);


% PCA for social
colorSoc = [repmat([0 0 0],[24,1]); repmat([1 0 0],[32,1]); repmat([.7 .7 .7],[20,1]); repmat([219 125 211]./256,[20,1])];

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca([BBsl; WWsl; BWsl; WBsl]); % 42.2598, 17.8366, 11.0766
figure(1); scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),200,colorSoc,'filled'); hold on; 
axis equal; 

% PCA for SOC
[COEFFS, SCORES, LATENTS, TSQUAREDS, EXPLAINEDS, MUS] = pca([BBs; WWs; BWs; WBs]); % 30.1104, 12.5784, 7.5877
figure(2); scatter3(SCORES(:,1),SCORES(:,2),SCORES(:,3),200,colorSoc,'filled'); hold on; 
axis equal; 

meanBBsl = mean(BBsl);
meanWWsl = mean(WWsl);
meanBWsl = mean(BWsl);
meanWBsl = mean(WBsl);

meanBBs = mean(BBs);
meanWWs = mean(WWs);
meanBWs = mean(BWs); 
meanWBs = mean(WBs);


figure(1); imagesc(log([BBsl; ones(2,118); WWsl; ones(2,118); BWsl; ones(2,118); WBsl])); colormap(gray); caxis([-9 -2]);
figure(2); imagesc(log([BBs; ones(2,133); WWs; ones(2,133); BWs; ones(2,133); WBs])); colormap(gray); caxis([-9 -2]);


figure(1); imagesc(log([BBsl; ones(2,118); BWsl; ones(2,118); WWsl; ones(2,118); WBsl])); colormap(gray); caxis([-9 -2]); axis off;
set(gcf,'Position',[278 815 669 172]);

figure(2); imagesc(log([BBs; ones(2,133); BWs; ones(2,133); WWs; ones(2,133); WBs])); colormap(gray); caxis([-9 -2]); axis off
set(gcf,'Position',[278 815 669 172]);



mBB = BBsl-meanBBsl; figure(1); imagesc(mBB); caxis([-.25 .25]); colormap(cmapdiff);
mWW = WWsl-meanBBsl; figure(2); imagesc(mWW); caxis([-.25 .25]); colormap(cmapdiff);
mWB = WBsl-meanBBsl; figure(3); imagesc(mWB); caxis([-.25 .25]); colormap(cmapdiff);
mBW = BWsl-meanBBsl; figure(4); imagesc(mBW); caxis([-.25 .25]); colormap(cmapdiff);

figure(1); imagesc([mBB; zeros(2,118); mBW; zeros(2,118); mWW; zeros(2,118); mWB]); colormap(cmapdiff); caxis([-.1 .1]);
set(gcf,'Position',[278 815 669 172]); axis off

figure(2); imagesc([mBBs; zeros(2,133); mBWs; zeros(2,133); mWWs; zeros(2,133); mWBs]); colormap(cmapdiff); caxis([-.1 .1]);
set(gcf,'Position',[278 815 669 172]); axis off




figure(1)
imagesc(log([meanBBsl; meanWWsl; meanBWsl; meanWBsl]))
axis off; colormap(gray); caxis([-9 -2.5]);
set(gcf,'Position',[210 885 696 91]);

figure(2)
imagesc(log([meanBBs; meanWWs; meanBWs; meanWBs]))
axis off; colormap(gray); caxis([-9 -2.5]);
set(gcf,'Position',[210 885 696 91]);


meanAll = (meanBBsl+meanWWsl+meanBWsl+meanBWsl)./4;
meanAllS = (meanBBs+meanWWs+meanBWs+meanBWs)./4;


mBB = BBsl-meanBBsl; figure(1); imagesc(mBB); caxis([-.25 .25]); colormap(cmapdiff);
mWW = WWsl-meanBBsl; figure(2); imagesc(mWW); caxis([-.25 .25]); colormap(cmapdiff);
mWB = WBsl-meanBBsl; figure(3); imagesc(mWB); caxis([-.25 .25]); colormap(cmapdiff);
mBW = BWsl-meanBBsl; figure(4); imagesc(mBW); caxis([-.25 .25]); colormap(cmapdiff);






