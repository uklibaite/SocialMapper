%% Code for sampling single movie for training set and re-embedding

% loading tracks, preprocessing
load('sampleTracksSocial.mat','p1','p2');
load('vecsVals_2022_12_22.mat','vecs','vals','L','C','mus');
pred1 = p1; pred2 = p2;
% minimal median filtering and smoothing
for x = 1:23
    for y = 1:3
        pred1(:,y,x) = smooth(medfilt1(p1(:,y,x),3),3);
        pred2(:,y,x) = smooth(medfilt1(p2(:,y,x),3),3);
    end
end

xIdx = 1:23; yIdx = 1:23;
[Xi Yi] = meshgrid(xIdx,yIdx);
Xi = Xi(:); Yi = Yi(:);
IDX = find(Xi~=Yi);
x1 = 4; x2 = 6; mf = 10; smf = 25;
nx = length(xIdx);
firstBatch = true;
currentImage = 0;
batchSize = 90000;

%% online PCA example
% given a cell 'predLONE' with many 23-keypoint trajectories (from lone and
% social context animals)
for j = 1:length(predLONE)
    try
        fprintf(1,['Processing batch # ' num2str(j) '\n']);
        ma1 = predLONE{j};
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
        
        scaleVal = lengtht(j)./250;
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
for j = 1:length(predSOC(:))
    try
        fprintf(1,['Processing batch # ' num2str(j) '\n']);
        ma1 = predSOC{j};
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
        
        scaleVal = lengthtS(j)./250;
        p1Dist = p1Dist.*scaleVal;
        
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
    catch
    end
end


L = currentImage; mu = mu ./ L; C = C ./ L - mu'*mu;
fprintf(1,'Finding Principal Components\n');
[vecs,vals] = eig(C); vals = flipud(diag(vals)); vecs = fliplr(vecs);
mus = mu;

% vecs/vals used for individual embedding
load('vecsVals_2022_12_22.mat','vecs','vals','L','C','mus');



vecs10 = vecs(:,1:10);
vecs15 = vecs(:,1:15);
minF = .5; maxF = 20; pcaModes = 20; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = 250;


ma1 = pred1;
sj = returnDist3d(squeeze(pred1(:,:,1)),squeeze(pred1(:,:,7)));
lz = prctile(sj,95);
scaleVal = 250./lz;

xIdx = 1:23; yIdx = 1:23;
[Xi Yi] = meshgrid(xIdx,yIdx);
Xi = Xi(:); Yi = Yi(:);
IDX = find(Xi~=Yi);
x1 = 4; x2 = 6; mf = 10; smf = 25; 
nx = length(xIdx);
firstBatch = true;
currentImage = 0;
batchSize = 90000;

load('train1228.mat')
trainingSetData = cD; trainingEmbeddingZ = ydata;

vecs15 = vecs(:,1:15);
minF = .5; maxF = 20; pcaModes = 20; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
fprintf(1,['training data = ' num2str(size(trainingSetData)) '\n']);
parameters.batchSize = 5000;


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
scaleVal = 250./lz;
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
p1z = (p1z+floorval).*scaleVal; 

nnData = [data2 .1*p1z' .5*jv];


fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(nnData,trainingSetData,...
    trainingEmbeddingZ,[],parameters);

z = zValues; z(~inConvHull,:) = zGuesses(~inConvHull,:);


