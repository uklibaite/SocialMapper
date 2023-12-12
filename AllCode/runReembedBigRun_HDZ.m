function [] = runReembedBigRun_HDZ(fileName,tData,vData,savePath)

addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

d1 = load(fileName);
pred = d1.pred;

load(tData);
load(vData);

trainingSetData = cD;
trainingEmbeddingCC = CC;
trainingEmbeddingZ = ydata;
K = 1;

pred = d1.pred;
for x = 1:23
    for y = 1:3
        pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
    end
end
ma1 = pred;
sj = returnDist3d(squeeze(pred(:,:,1)),squeeze(pred(:,:,7)));
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
batchSize = 5000;
N = length(nnData(:,1));
batches = ceil(N/batchSize);
readout = 2000;
cGuesses = zeros(N,1);
 
for j = 1:batches
tic
fprintf(1,'\t Processing batch #%4i out of %4i\n',j,batches);
idx = (1:batchSize) + (j-1)*batchSize;
idx = idx(idx<=N);
current_guesses = zeros(length(idx),1);
currentData = nnData(idx,:);
D2 = pdist2(currentData,trainingSetData,'cityblock');

[Dt Dt2] = sort(D2,2,'ascend');
Dt = Dt(:,1:K); Dt2 = Dt2(:,1:K); DtC = trainingEmbeddingCC(Dt2);

mddt = mode(DtC,2);
current_guesses(:,1) = mddt;
cGuesses(idx,:) = current_guesses;
toc
end


tic
fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(nnData,trainingSetData,...
    trainingEmbeddingZ,[],parameters);

z = zValues; z(~inConvHull,:) = zGuesses(~inConvHull,:);
toc

ss = strsplit(fileName,'/');
spath = [savePath ss{end-4} '_' ss{end-3} '_' ss{end-1}(end-3:end) '_RE.mat'];
save(spath,'z','inConvHull','cGuesses','fileName','parameters');