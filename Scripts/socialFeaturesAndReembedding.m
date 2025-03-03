%% Code for sampling single movie for training set and re-embedding

% loading tracks, preprocessing
load('sampleTracksSocial.mat','p1','p2');
load('vecsValsSOC_2022_12_28.mat','vecs','vals','L','C','mus');
pred1 = p1; pred2 = p2;
% minimal median filtering and smoothing
for x = 1:23
    for y = 1:3
        pred1(:,y,x) = smooth(medfilt1(p1(:,y,x),3),3);
        pred2(:,y,x) = smooth(medfilt1(p2(:,y,x),3),3);
    end
end

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

vecs6 = vecs(:,1:6);
minF = .2; maxF = 5; pcaModes = 10; numModes = pcaModes;
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.minF = minF;
parameters.maxF = maxF;
numPerDataSet = 250;

% perspective of animal 1
ma1 = pred1;
ma2 = pred2;
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
nnData = [data2 .1*rat1z .1*rat2z p1v(:,[1 11 12 15 19 23]) p2v(:,[1 11 12 15 19 23]) 10*angInfo .05*PJ 2*PJN];
yData = tsne(nnData(1:20:end,:));
amps = sum(nnData',1);
dataMat1 = nnData;

% FIND GOOD TEMPLATES AND SAVE THEM
[signalData,signalAmps] = findTemplatesFromData(...
    nnData(1:20:end,:),yData,amps(1:20:end,:),numPerDataSet,parameters);
% mDC{rec,1} = signalData; mDCA{rec,1} = signalAmps;

% perspective of animal 2
ma1 = pred2;
ma2 = pred1;
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

nnData = [data2 .1*rat2z .1*rat1z p2v(:,[1 11 12 15 19 23]) p1v(:,[1 11 12 15 19 23]) 10*angInfo2 .05*PJ 2*PJN];
yData = tsne(nnData(1:20:end,:));
amps = sum(nnData',1);
dataMat2 = nnData;
%         % FIND GOOD TEMPLATES AND SAVE THEM
[signalData,signalAmps] = findTemplatesFromData(...
    nnData(1:20:end,:),yData,amps(1:20:end,:),numPerDataSet,parameters);
% mDC{rec,2} = signalData; mDCA{rec,2} = signalAmps;

%fileName = FNSX{rec};
%ss = strsplit(fileName,'/');
% save(['socialSaves/' ss{end-4} 'x' ss{end-3} '.mat'],'dataMat1','dataMat2','fileName');

%% Reembedding into joint social map

load('trainingSocial_20230131.mat')
trainingSetData = trainingData;

parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.batchSize = 5000;

fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(dataMat1,trainingSetData,...
    trainingEmbedding,[],parameters);

z1 = zValues; z1(~inConvHull,:) = zGuesses(~inConvHull,:);
inCH1 = inConvHull;

fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(dataMat2,trainingSetData,...
    trainingEmbedding,[],parameters);

z2 = zValues; z2(~inConvHull,:) = zGuesses(~inConvHull,:);
inCH2 = inConvHull;

% z1 and z2 are 2-d trajectories through low-dimensional joint
% behavior space
