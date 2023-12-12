function [] = runReembedBigRunS_CCZ(fileName,tData,savePath)

addpath(genpath('./utilities/'));
addpath(genpath('./tSNE/'));

fname1 = fileName;
load(fileName);
load(tData);

trainingSetData = trainingData;
K = 1; 
parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.batchSize = 5000;

fprintf(1,'Finding Embeddings\n');
batchSize = 5000;
N = length(dataMat1(:,1));
batches = ceil(N/batchSize);
readout = 2000;
cGuessesRat1 = zeros(N,1);

for j = 1:batches
    tic
    fprintf(1,'\t Processing batch #%4i out of %4i\n',j,batches);
    idx = (1:batchSize) + (j-1)*batchSize;
    idx = idx(idx<=N);
    current_guesses = zeros(length(idx),1);
    currentData = dataMat1(idx,:);
    D2 = pdist2(currentData,trainingSetData,'euclidean');
    
    [Dt Dt2] = sort(D2,2,'ascend');
    Dt = Dt(:,1:K); Dt2 = Dt2(:,1:K); DtC = CC(Dt2);
    
    mddt = mode(DtC,2);
    current_guesses(:,1) = mddt;
    cGuessesRat1(idx,:) = current_guesses;
    toc
end

fprintf(1,'Finding Embeddings\n');
batchSize = 5000;
N = length(dataMat2(:,1));
batches = ceil(N/batchSize);
readout = 2000;
cGuessesRat2 = zeros(N,1);

for j = 1:batches
    tic
    fprintf(1,'\t Processing batch #%4i out of %4i\n',j,batches);
    idx = (1:batchSize) + (j-1)*batchSize;
    idx = idx(idx<=N);
    current_guesses = zeros(length(idx),1);
    currentData = dataMat2(idx,:);
    D2 = pdist2(currentData,trainingSetData,'euclidean');
    
    [Dt Dt2] = sort(D2,2,'ascend');
    Dt = Dt(:,1:K); Dt2 = Dt2(:,1:K); DtC = CC(Dt2);
    
    mddt = mode(DtC,2);
    current_guesses(:,1) = mddt;
    cGuessesRat2(idx,:) = current_guesses;
    toc
end



parameters = setRunParameters([]);
parameters.samplingFreq = 50;
parameters.batchSize = 5000;

tic
fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(dataMat1,trainingSetData,...
    trainingEmbedding,[],parameters);

z1 = zValues; z1(~inConvHull,:) = zGuesses(~inConvHull,:);
inCH1 = inConvHull;
toc

tic
fprintf(1,'Finding Embeddings\n');
[zValues,zCosts,zGuesses,inConvHull,meanMax,exitFlags] = ...
    findTDistributedProjections_fmin(dataMat2,trainingSetData,...
    trainingEmbedding,[],parameters);

z2 = zValues; z2(~inConvHull,:) = zGuesses(~inConvHull,:);
inCH2 = inConvHull;
toc


ss = strsplit(fname1,'/');
spath = [savePath ss{end}(1:end-4) '_SRE.mat'];
save(spath,'z1','z2','inCH1','inCH2','cGuessesRat1','cGuessesRat2','fileName','parameters');