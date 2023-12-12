%% 
addpath(genpath('/home/ugne/Dropbox/MultiFlyAnalysis_OLD/tSNE'))
skeleton = load('/home/ugne/Label3D/skeletons/rat23.mat')

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

skeleton.color = scM;%zeros(23,3); 
skeleton.joints_idx = joints_idx;
skeleton.mcolor = scM;
% skeleton.joint_names = joint_names;



%% 





load('A_SDANNCE/bigRun/results202302112.mat','isSOC','isAMPH','ratGROUP',...
    'ratID','ratGID','ratDATE','ratCZZ','ratCZ','ratPID','ratPGID',...
    'ratSOCEZ','ratSOCEZZ','ratCLASSHIST_Z','coarseCLASSHIST_Z',...
    'ratSOCHIST_Z','coarseSOCHIST_Z','ratTC','ratTS','ratTCself',...
    'ratTSself');

load('A_SDANNCE/bigRun/results202302112info.mat','socOrderZ','socOrderTagZ',...
    'socColors','orderedCZ','orderedCZCoarse','gS','gSA','gSLE','FILECOMBO',...
    'mdistC','ratGEN','ratPGEN');

load('/home/ugne/Dropbox/A_SDANNCE/bigRun/classicalEmbeddingSaves.mat', 'LL','LL2')
load('/home/ugne/Dropbox/A_SDANNCE/bigRun/classicEmbeddingInfo.mat','colorsCoarse','llbwb','D','LLLC','labelsCZCoarse','labelsCZ');


% individual 
xx = linspace(-58,58,501);
vSmooth = .5;
medianLength = 1;
pThreshold = [];
minRest = 5;
obj = [];
fitOnly = true;
numGMM = 2;

evall = combineCells(ratCZZ);
[wr,segments,v,obj,pRest,vals,vx,vy] = ...
    findWatershedRegions_v2(evall,xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);

wrALL = cell(size(ratCZZ));
for i = 1:1002
    try
    [wrALL{i},segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(ratCZZ{i},xx,LL,vSmooth,medianLength,pThreshold,minRest,obj,false,numGMM);
    catch
    end
end

%labelsCZCoarse = behLabels; 
labelsCZCoarse(61) = 1;
wFINE = cell(size(wrALL)); wCOARSE = cell(size(wrALL));
for i = 1:length(wrALL)
    w = wrALL{i};
    w(find(w==0)) = nan;
    w = fillmissing(w,'nearest');
    wFINE{i} = w; 
    
    w = wrALL{i};
    w(w~=0) = labelsCZCoarse(w(w~=0));
    w(find(w==0)) = nan;
    w = fillmissing(w,'nearest');
    wCOARSE{i} = w;
end


% final sorting for plots
allCV = cell(1002,1);
for i = 1:1002
    try
        pred = predALL{i,1};
        jd= [0; smooth(sqrt(sum(diff(squeeze(pred(:,:,7))).^2,2)),20)];
        allCV{i} = jd;
    catch
    end
end
allcv = combineCells(allCV);
allev = combineCells(wrALL);

meancv = zeros(1,163);
for b = 1:163
    idx = find(allev==b);
    mdc = median(allcv(idx),'omitnan');
    meancv(b) = mdc;
end


% meanpower = zeros(size(meancv));
% for i = 1:163
%     meanpower(i) = sum(sum(W_xyz(:,:,i)));
% end


[groups,~,~] = makeGroupsAndSegments(wrALL(:),max(max(LL)),ones(1,1002),20);

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
    
load('updatedSaves/individualClusteringInfo.mat','behOrderTagZ','behOrderZ','behLabels')

imagesc(LL2); hold on;
for i = 1:163
    [xi yi] = find(LL2==i);
    mxi = mean(xi); myi = mean(yi);
    text(myi,mxi,num2str(i));
end

LL3 = zeros(size(LL2));
for i = 1:163
    cid = find(LL2==i);
    LL3(cid) = behLabels(i);
end
imagesc(LL3); colormap(colorsCoarse);

fid = fopen('bigRun/b_CLASS_Z_labels.txt','r');
fSZ = textscan(fid, "%s", "Delimiter", "\n");
fclose(fid);
bdef = fSZ{1};
bdef3 = cell(163,3);
for i = 1:163
    bd = strsplit(bdef{i},' ');
    bdef3{i,1} = bd{1};
    bdef3{i,2} = bd{2};
    bdef3{i,3} = bd{3};
end
behLabels = nan(163,1);
for i = 1:163
    behLabels(i) = str2num(bdef3{i,2});
end

behOrderZ = []; behOrderTagZ = [];
for b = 1:9
    xi = find(behLabels==b);
    dxi = meancv(xi);
    [dsorted,dsxi] = sort(dxi,'ascend');
    sortedxi = xi(dsxi);
    behOrderZ = [behOrderZ; sortedxi];
    behOrderTagZ = [behOrderTagZ; b*ones(size(sortedxi))];
end

b = bar(meancv(behOrderZ));
b.FaceColor = 'flat';
for i = 1:length(behOrderZ)
    b.CData(i,:) = colorsCoarse(behOrderTagZ(i)+1,:);
end

% individual: behOrderTagZ, behOrderZ, behLabels || social: socOrderTagZ,
% socOrderZ, lsocZtag

socOrder = [0 3 1 4 2 7 5 6 8];
socOrderZ = []; socOrderTagZ = [];
for i = 1:9
    ti = socOrder(i);
    xi = find(lsocZtag2==ti);
    dxi = meandc2(xi);
    [sdsorted,sdxi] = sort(dxi,'descend');
    sortedxi = xi(sdxi);
    socOrderZ = [socOrderZ; sortedxi];
    socOrderTagZ = [socOrderTagZ; ti*ones(size(sortedxi))];
end

socColors2 = socColors; socColors2(2,:) = [.8 .8 .8];

b = bar(meandc2(socOrderZ));
b.FaceColor = 'flat';
for i = 1:length(socOrderZ)
    b.CData(i,:) = socColors2(socOrderTagZ(i)+2,:);
end
saveas(gcf,'figs20230910/barsMeanDC_sb.tif')

axis off 
saveas(gcf,'figs20230910/barsMeanDC.tif')

figure(2);
b2 = bar(shwr);
b2.FaceColor = 'flat';
for i = 1:socN
    bid = shwridx(i);
    id = find(socOrderZ==bid);
    ccol = socColors2(socOrderTagZ(id)+2,:);
    b2.CData(i,:) = ccol;
end

[socOrderTagZ socOrderZ]

load('updatedSaves/socialClusterInfo.mat','socOrderTagZ','socOrderZ','lsocZtag','lsocZtag2','socColors2','meandc','meandc2')

x0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0);
x0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1);
x1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0);
x1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1);
x2l = find(ratGROUP==6 & isAMPH==2 & isSOC==0);
x2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1);

zx0l = combineCells(ratCZZ(x0l));
zx0s = combineCells(ratCZZ(x0s));
[dx,dx0l] = findPointDensity(zx0l,1,501,[-58 58]);
[dx,dx0s] = findPointDensity(zx0s,1,501,[-58 58]);


[D,dx] = findPointDensity(combineCells(ratCZZ(:,1)),1,501,[-58 58]);


parameters = setRunParameters([]);
parameters.pcaModes = 69;
parameters.samplingFreq = 50;
parameters.minF = .25;
parameters.maxF = 15;
numModes = 69;
numW = 25;
BIG_W = cell(60,1);
L_W = cell(60,1);
maxNum = 163; N = 163;
minLength = 10;

POWER = zeros(163,575); powerN = zeros(163,1);

for rec = 1:1002
    rec
    try
        
        p1 = predALL{rec,1};
        p2 = reshape(p1,[90000 69]);
        [W1 a] = findWavelets(p2,numModes,parameters);
        wrt = wrALL{rec};
        G1 = makeGroupsAndSegments({wrt},maxNum,1,minLength);
        
        
 BIdx = cell(1,N);
    SortedW1 = cell(size(BIdx));
    for i = 1:N
        gtemp = G1{i};
        btemp = [];
        SW1t = cell(size(gtemp,1),1);
        for j = 1:size(gtemp,1)
            btemp = [btemp gtemp(j,2):gtemp(j,3)];
            SW1t{j} = W1(gtemp(j,2):gtemp(j,3),:);
        end
        BIdx{i} = btemp;
        SortedW1{i} = SW1t;
    end
    
    allW = cell(size(BIdx));
    for i = 1:N
        allW{i} = W1(BIdx{i},:);
    end
    
    allW = cell(size(BIdx));
    allP = cell(size(BIdx));
    for i = 1:N
        allW{i} = W1(BIdx{i},:);
        allP{i} = p2(BIdx{i},:);
    end
    
    
    allWM = zeros(69,25,N);
    l_w = zeros(1,N);
    for i = 1:N
        tempwm = mean(allW{i},1);
        l_w(i) = size(allW{i},1);
        allWM(:,:,i) = reshape(tempwm,[25 69])';
    end

    L_W{rec} = l_w;
    BIG_W{rec} = allWM;
    catch
    end
end

W_all = zeros(69,25,N);
L_all = zeros(1,N);
for i = 1:N
    TW = zeros(69,25);
    TL = 0;
    for M = 1:133
        try
        tw = BIG_W{M}(:,:,i);
        tl = L_W{M}(i);
        if sum(sum(isnan(tw)))==0
            TW = TW+tw;
            TL = TL+tl;
        end
        catch
        end
    end
    W_all(:,:,i) = TW./TL;
    L_all(i) = TL;
end

for i= 1:N
    WN(:,:,i) = W_all(:,:,i)./sum(sum(W_all(:,:,i)));
end

W_xyz = zeros(23,25,N);
for i = 1:N
    tw = W_all(:,:,i);
    tw2 = tw(1:3:end,:)+tw(2:3:end,:)+tw(3:3:end,:);
    W_xyz(:,:,i) = tw2;
end

WN_xyz = zeros(23,25,N);
for i = 1:N
    tw = WN(:,:,i);
    tw2 = tw(1:3:end,:)+tw(2:3:end,:)+tw(3:3:end,:);
    WN_xyz(:,:,i) = tw2;
end

% social
load('socialEmbeddingInfo2023_03_14.mat','FNSX','LL2soc','LLsoc','d','llbwbsoc','socColors')
load('updatedSaves/socialClusterInfo.mat','socOrderTagZ','socOrderZ','lsocZtag','lsocZtag2','socColors2','meandc','meandc2')

xx = linspace(-75,75,501);
vSmooth = .5;
medianLength = 1;
pThreshold = [];
minRest = 10;
obj = [];
fitOnly = true;
numGMM = 2;

evall = combineCells(ratSOCEZZ);
[wr,segments,v,obj,pRest,vals,vx,vy] = ...
    findWatershedRegions_v2(evall,xx,LLsoc,vSmooth,medianLength,pThreshold,minRest,obj,fitOnly,numGMM);

wrsocALL = cell(size(ratSOCEZZ));
for i = 1:1002
    try
    [wrsocALL{i},segments,v,obj,pRest,vals,vx,vy] = ...
            findWatershedRegions_v2(ratSOCEZZ{i},xx,LLsoc,vSmooth,medianLength,pThreshold,minRest,obj,false,numGMM);
    catch
    end
end

wsocFINE = cell(size(ratCZ)); wsocCOARSE = cell(size(ratCZ));
for i = 1:length(wrsocALL)
    w = wrsocALL{i};
    w(find(w==0)) = nan;
    w = fillmissing(w,'nearest');
    wsocFINE{i} = w; 
    
    w(w~=0) = lsocZtag(w(w~=0));
    wsocCOARSE{i} = w;
end


save('updatedSaves/allClustersFineCoarse.mat','wFINE','wCOARSE','wsocFINE','wsocCOARSE');


socialConstruct = [1 1 1 1 1 1 1 1 1; 
    1 1 4 1 4 4 4 4 4;
    1 5 1 5 8 5 8 8 8;
    1 1 4 1 2 2 2 2 2;
    1 5 8 3 1 3 3 3 3;
    1 5 4 3 2 8 8 8 8;
    1 5 8 3 2 8 8 8 8;
    1 5 8 3 2 8 8 8 9;
    1 5 8 3 2 8 8 9 9];


socialMatrix = nan(161,161);
for b = 1:161
    for j = 1:161
        cons = socialConstruct(lsocZtag(b)+1,lsocZtag(j)+1);
        socialMatrix(b,j) = cons;
    end
end
SMAT = nan(size(socialMatrix));
SMAT(socialMatrix==1) = 1;
SMAT(socialMatrix==2) = 2;
SMAT(socialMatrix==3) = 3;
SMAT(socialMatrix==4) = 2;
SMAT(socialMatrix==5) = 3;
SMAT(socialMatrix==6) = 4;
SMAT(socialMatrix==7) = 4;
SMAT(socialMatrix==8) = 4;
SMAT(socialMatrix==9) = 4;


wsocSYMM = cell(size(ratCZ));
SMSOC = cell(1002,1);

for i = 1:396
    w1 = wsocFINE{210+i};
    w2 = wsocFINE{210+i+396};
    smat1 = zeros(size(w1)); smat2 = zeros(size(w2));
    for n = 1:length(w1)
        smat1(n) = SMAT(w1(n),w2(n)); smat2(n) = SMAT(w2(n),w1(n));
    end
    SMSOC{210+i} = smat1;
    SMSOC{210+i+396} = smat2;
end





% LONGEVANS_M_SOC6/2022_11_28_M1_M4/, REC=388
zh1 = squeeze(predALL{210+178}(53001:58000,3,[1 3 4 5 7 11 15 19 23]));
zh2 = squeeze(predALL{210+178+396}(53001:58000,3,[1 3 4 5 7 11 15 19 23]));

figure(1);
for i = 1:9
    plot(1:5000,zh1(:,i)','Color','k','LineWidth',1); hold on;
    plot(1:5000,zh2(:,i)'-350,'Color','k','LineWidth',1); hold on;
end

set(gcf,'Position',[75 674 1846 346]); axis off
saveas(gcf,'figs20230910/combo_z.tif');
saveas(gcf,'figs20230910/combo_z.eps');

cw_wt = zeros(55,90000);
t  = 178;
wt1 = wCOARSE{210+t}; wt2 = wCOARSE{210+t+396};
bid1 = 2:6:55; 
bid2 = 4:6:55;
for b = 1:9
    bi = find(wt1==b);
    cw_wt(bid1(b),bi) = b;
    bi = find(wt2==b);
    cw_wt(bid2(b),bi) = b;
end

imagesc(cw_wt); colormap(colorsCoarse); axis ij
axis([53000 58000 0 55]);



% SCN2A_SOC1/2022_09_22_M1_M2/, REC = 342
zh1 = squeeze(predALL{342}(5001:8000,3,[1 3 4 5 7 11 15 19 23]));
zh2 = squeeze(predALL{342+396}(5001:8000,3,[1 3 4 5 7 11 15 19 23]));

figure(1);
for i = 1:9
    plot(1:3000,zh1(:,i)','Color','k','LineWidth',1); hold on;
    plot(1:3000,zh2(:,i)'-350,'Color','k','LineWidth',1); hold on;
end






set(gcf,'Position',[75 674 1846 346]); axis off
saveas(gcf,'figs20230910/combo_z_touch.tif');
saveas(gcf,'figs20230910/combo_ztouch.eps');


REC = 342;
WR_social1 = ratSOCEZ{REC};
WR_social2 = ratSOCEZ{REC+396};

f = touchFiles{REC-210};
contacts = h5read(f,'/contacts');
fnums = h5read(f,'/frames');

contactsV = double(contacts+1);
contactsT = vertex_idxRed(contactsV);
zeroidx1 = find(contactsT(1,:)==0);
zeroidx2 = find(contactsT(2,:)==0);
zeroidx = [zeroidx1 zeroidx2];

contactsT(:,zeroidx) = [];
contactsV(:,zeroidx) = [];
fnumsF = fnums+1;
fnumsF(zeroidx) = [];

h1 = hist(contactsV(1,:),1:6880);
h2 = hist(contactsV(2,:),1:6880);

n = 90000;
tactogram = zeros(49,n);
tic
for i = 1:n
    fx = find(fnumsF==i);
    allcred = contactsT(:,fx)';
    ucred = unique(allcred,'rows');
    temptouch = zeros(7,7);
    temptouch(ucred(:,1),ucred(:,2)) = 1;
    tactogram(:,i) = temptouch(:);
end
toc

tactIdx = tactogram(:,5001:8000);
tactIdx2 = zeros(size(tactIdx));
for i = 1:49
    t1 = tactIdx(i,:);
    t2 = round(medfilt1(tactIdx(i,:),5));
    cc = bwconncomp(t2);
    for j = 1:length(cc.PixelIdxList)
        if length(cc.PixelIdxList{j})>9
                tactIdx2(i,cc.PixelIdxList{j}) = 1;
        end
    end
    
end

imagesc(tactIdx+tactIdx2);
figure(2); imagesc(-tactIdx2); colormap(gray); axis off
saveas(gcf,'figs20230910/tactogram3001ax.tif')
colSocialReduced = [.9 .9 .9; 0 0 1; 1 0 0; .6 0 .6];
r1c = wCOARSE{REC};
r2c = wCOARSE{REC+396};
swr = SMSOC{REC};

figure(1); imagesc(r1c(5001:8000)'); set(gcf,'Position',[1 974 1891 56]); colormap(colorsCoarse); caxis([0 9]); axis off; %saveas(gcf,'figs20230910/rat1c_5001.tif');
figure(2); imagesc(r2c(5001:8000)'); set(gcf,'Position',[1 974 1891 56]); colormap(colorsCoarse); caxis([0 9]); axis off; %saveas(gcf,'figs20230910/rat2c_5001.tif');
figure(3); imagesc(swr(5001:8000)'); set(gcf,'Position',[1 974 1891 56]); colormap(colSocialReduced); caxis([1 4]); axis off;  %saveas(gcf,'figs20230910/ratSOC_5001.tif');

%% building data struct and saving code
grouplist = [{'ARID1B'},{'CHD8'},{'GRINB'},{'NRXN1'},{'SCN2A'},{'LONGEVANS'},{'ratsInColor'}];
gidscale = [100; 200; 300; 400; 500; 600; 700];

arid1bG = [0 1 1 1 0 0 1 1 0 0 1 0];
chd8G = [0 0 1 1 0 0 1 1];
grinbG = [1 1 0 1 0 0];
nrxn1G = [0 0 1 1 0 1 0 1];
scn2aG = [0 0 0 1 1 1];
cntnapG = [.5 1 1 2 0 0 0]; % wt/het/hom/female

arid1bCAGE_G = [1 0 0 0 1 1 0 0 1 1 0 1];
chd8CAGE_G = [0 0 1 1 0 0 1 1];
grinbCAGE_G = [0 0 1 0 1 1];
nrxn1CAGE_G = [0 0 1 1 1 0 1 0];
scn2aCAGE_G = [1 1 1 0 0 0];

arid1bLookup = 100 + (1:12);
chd8Lookup = 200 + (1:8);
grinbLookup = 300 + (1:6);
nrxn1Lookup = 400 + (1:8);
scn2aLookup = 500 + (1:6);
leLookup = 600 + (0:5);
cntnapLookup = 700 + (1:7);


amph_le_soc6 = [0 0; 0 0; 0 0; 1 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 1 0; 0 1; 0 1; 0 0; 1 0];
amph_le_soc7 = [0 1; 1 0; 0 1; 0 0; 0 0; 0 0; 1 0; 1 0; 0 0; 0 1; 0 1; 0 0; 1 0; 0 0; 0 0];
amph_le_soc8 = [0 1; 1 0; 0 1; 0 0; 0 0; 0 0; 0 1; 0 1; 0 0; 1 0; 0 0; 0 0; 0 0; 0 1; 1 0];
amphLE = [amph_le_soc6; amph_le_soc7; amph_le_soc8];
amphLE2 = amphLE;
for i = 1:length(amphLE)
    if amphLE(i,1)==1
        amphLE(i,2)=2;
    else
    end
    if amphLE(i,2)==1
        amphLE(i,1)=2;
    else
    end
end


% FILL IN ALL amphLE - NEED LE_WEEK_8!! 


% 210 LONE, 396x2 SOCIAL
% 
isamphL = zeros(210,1); isamphS = zeros(396,2);
isamphL(119:124) = 1;
isamphS(177:221,:) = amphLE;









% save kinematics, lone embedding, social embedding, .. all ratID-style
% info

% include filename (extension only)
% 


FILECOMBO = [FN; FNS(:,1); FNS(:,2)];

load('A_SDANNCE/bigRun/fileNames_1002.mat','FNX_lone','FNX_soc');


FN_lone = cell(size(FNX_lone));
FN_soc = cell(size(FNX_soc));

for i = 1:210
    FN_lone{i} = FNX_lone{i}(111:end);
end

for i = 1:396
    FN_soc{i,1} = FNX_soc{i,1}(111:end);
    FN_soc{i,2} = FNX_soc{i,2}(111:end);
end

projectFolders = cell(size(FILECOMBO));
for i = 1:length(FILECOMBO)
    t = FILECOMBO{i}(111:end);
    t2 = strsplit(t,'/');
    projectFolders{i} = [t2{1} '/' t2{2}];
    movieF{i} = [base t2{1} '/' t2{2} '/videos/Camera1/0.mp4'];
end



colorWT = [131 155 242
232 200 121
191 132 209
84 143 84
230 127 119];
colorKO = [12 62 242
252 184 10
174 24 219
17 153 17
242 36 22];



base = '/run/user/1000/gvfs/smb-share:server=holy-isilon.rc.fas.harvard.edu,share=olveczky_lab/Lab/dannce_rig/';

save([base 'test.mat'],pred);



for i = 1:396
    i
    fn = []; pred = []; clear d1
    try
        fn1 = [base FN_soc{i,1}]; d1 = load(fn1); pred1 = d1.pred;
        fn2 = [base FN_soc{i,2}]; d2 = load(fn2); pred2 = d2.pred;
        
        sdannce_struct.rat1.pred = pred1;
        sdannce_struct.rat2.pred = pred2;
        
        sdannce_struct.rat1.ratGROUP = grouplist{ratGROUP(i+210)};
        sdannce_struct.rat2.ratGROUP = grouplist{ratGROUP(i+210)};
        
        sdannce_struct.rat1.ratID = ratID{i+210};
        sdannce_struct.rat2.ratID = ratPID{i+210};
        
        sdannce_struct.rat1.isSOC = isSOC(i+210);
        sdannce_struct.rat2.isSOC = isSOC(i+210+396);
        
        sdannce_struct.rat1.isAMPH = isAMPH(i+210);
        sdannce_struct.rat2.isAMPH = isAMPH(i+210+396);
        
        sdannce_struct.rat1.ratDATE = ratDATE(i+210);
        sdannce_struct.rat2.ratDATE = ratDATE(i+210+396);
        
        sdannce_struct.rat1.autoWR = ratCZ{210+i};
        sdannce_struct.rat2.autoWR = ratCZ{210+i+396};
        
        sdannce_struct.rat1.autoWRC = ratCZ{210+i};
        sdannce_struct.rat2.autoWRC = ratCZ{210+i+396};
        
        sdannce_struct.rat1.socWR = ratSOCEZ{i+210};
        sdannce_struct.rat2.socWR = ratSOCEZ{i+210+396};
        
        sdannce_struct.rat1.ratGEN = ratGEN(i+210);
        sdannce_struct.rat2.ratGEN = ratGEN(i+210+396);
        
        sdannce_struct.rat1.ratPGEN = ratPGEN(i+210);
        sdannce_struct.rat2.ratPGEN = ratPGEN(i+210+396);
    catch
    end
end

%h5create('A_SDANNCE/test.h5',"/rat1",[20 Inf],"Chunksize",[5 5])
%h5write('A_SDANNCE/test.h5','/rat1',pred)

        %sdannce_struct.rat1.autoEZ = ratCZZ{i+210};
        %sdannce_struct.rat2.autoEZ = ratCZZ{i+210+396};

        %sdannce_struct.rat1.socEZ = ratSOCEZZ{i+210};
        %sdannce_struct.rat2.socEZ = ratSOCEZZ{i+210+396};
        
        
isSOC = [isSOCL; isSOCS(:)];
isAMPH = [isamphL; isamphS(:)];
ratGROUP = [ratGROUP_lone; ratGROUP_soc(:)];
ratID = [ratID_lone; ratID_soc(:)];
ratGID = [ratGID_lone; ratGID_soc(:)];
ratDATE = [ratDATE_lone; ratDATE_soc(:)];
ratCZZ = [ratCZZ_lone; ratCZZ_soc(:)];
ratCZ = [ratCZ_lone; ratCZ_soc(:)];
%ratCC = [ratCC_lone; ratCC_soc(:)];
ratPID = [ratPID_lone; ratPID_soc(:)];
ratPGID = [ratPGID_lone; ratPGID_soc(:)];
%ratSOCEC = [ratSOCEC_lone; ratSOCEC_soc(:)];
ratSOCEZ = [ratSOCEZ_lone; ratSOCEZ_soc(:)];
ratSOCEZZ = [ratSOCEZZ_lone; ratSOCEZZ_soc(:)];




% if from wt/wt or ko/ko group - save in separate folder
% 



load('/home/ugne/Dropbox/A_SDANNCE/bigRun/social_embedding_info_20230418.mat')


ls418 = lsocZtag;
load('/home/ugne/Dropbox/A_SDANNCE/socialEmbeddingInfo2023_03_14.mat')

ls314 = lsocZtag;

LL3soc = -1*ones(size(LL2soc));
for i = 1:161
    LL3soc(find(LL2soc==i)) = lsocZtag(i);
end

load('/home/ugne/Dropbox/A_SDANNCE/socialEmbeddingInfo2023_03_14.mat')
LL3soc = -1*ones(size(LL2soc));
for i = 1:161
    LL3soc(find(LL2soc==i)) = lsocZtag(i);
end


imagesc(LL3soc); axis xy; colormap(socColors2);
axis off equal
saveas(gcf,'figs20230910/LL3.tif')

for i = 1:161
    [xi yi] = find(LL2soc==i);
    mxi = mean(xi); myi = mean(yi);
    text(myi,mxi,num2str(i));
end
saveas(gcf,'figs20230910/LL3numbered.tif')

%% reading from shared data 
allDC = cell(396,1);
for i = 1:396
    try
    %fn1 = [base FN_soc{i,1}]; d1 = load(fn1); ma1 = d1.pred;
    %fn2 = [base FN_soc{i,2}]; d2 = load(fn2); ma2 = d2.pred;
    ma1 = predALL{210+i,1}; ma2 = predALL{210+i,2};
    p1_xy = ma1(:,:,[1 5 7]);
    p2_xy = ma2(:,:,[1 5 7]);
    n1 = p1_xy(:,:,1); c1 = p1_xy(:,:,2); t1 = p1_xy(:,:,3);
    n2 = p2_xy(:,:,1); c2 = p2_xy(:,:,2); t2 = p2_xy(:,:,3);
    distc = sqrt(sum((c1-c2).^2,2));
    allDC{i} = distc;
    catch
    end
end
alldistc = [combineCells(allDC); combineCells(allDC)];

allsocwr = combineCells(ratSOCEZ(211:end));

hwr = hist(allsocwr,0:161);
hwr2 = hwr(2:end);
[shwr shwridx] = sort(hwr2);

tfsoc = sum(hwr);

meandc = zeros(1,161);
for b = 1:161
    idx = find(allsocwr==b);
    mdc = mean(alldistc(idx),'omitnan');
    meandc(b) = mdc;
end

meandc2 = meandc;
meandc2(103) = NaN; meandc2(129) = NaN; meandc2(125) = NaN;
fid = fopen('bigRun/b_SOCL_SEPT9.txt','r');
fSZ = textscan(fid, "%s", "Delimiter", "\n");
fclose(fid);

lsocZ = fSZ{1};
socN = 161;
lsocZtag = zeros(socN,1);
for i = 1:socN
    tag = str2double(lsocZ{i});
    lsocZtag(i) = tag;
end

lsocZtag2 = lsocZtag; 
lsocZtag2(1) = -1; lsocZtag2(94) = -1; lsocZtag2(103) = -1; lsocZtag2(129) = -1; lsocZtag2(125) = -1;
% socOrder = [0 3 4 1 2 7 5 6 8];
socOrder = [0 3 1 4 2 7 5 6 8];
socOrderZ = []; socOrderTagZ = [];
for i = 1:9
    ti = socOrder(i);
    xi = find(lsocZtag2==ti);
    dxi = meandc2(xi);
    [sdsorted,sdxi] = sort(dxi,'descend');
    sortedxi = xi(sdxi);
    socOrderZ = [socOrderZ; sortedxi];
    socOrderTagZ = [socOrderTagZ; ti*ones(size(sortedxi))];
end

socColors2 = socColors; socColors2(2,:) = [.8 .8 .8];

b = bar(meandc2(socOrderZ));
b.FaceColor = 'flat';
for i = 1:length(socOrderZ)
    b.CData(i,:) = socColors2(socOrderTagZ(i)+2,:);
end
saveas(gcf,'figs20230910/barsMeanDC_sb.tif')

axis off 
saveas(gcf,'figs20230910/barsMeanDC.tif')

figure(2);
b2 = bar(shwr);
b2.FaceColor = 'flat';
for i = 1:socN
    bid = shwridx(i);
    id = find(socOrderZ==bid);
    ccol = socColors2(socOrderTagZ(id)+2,:);
    b2.CData(i,:) = ccol;
end

[socOrderTagZ socOrderZ]

save('updatedSaves/socialClusterInfo.mat','socOrderTagZ','socOrderZ','lsocZtag','lsocZtag2','socColors2','meandc','meandc2')

%%

[groups40,~,~]= makeGroupsAndSegments(ratSOCEZ,161,ones(1002,1),40);
[groups25,~,~]= makeGroupsAndSegments(ratSOCEZ,161,ones(1002,1),25);
lg = zeros(161);
for i = 1:161
    lg(i) = length(groups40{i});
end

offx = 800;
offy = 800;
offz = 0;
coffx = zeros(4,4); coffy = zeros(4,4); coffz = zeros(4,4);
for i = 1:4
    for j = 1:4
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy';



isrr = ones(1,1002);
for tm = 1:161
    try
        if meandc(tm)>300
        G = groups40{tm};
        lgs = size(G,1);
        if lgs>16
            newG = G(randperm(lgs,16),:);
        else
            newG = G;
        end
        lenG = 150;
        gMarkers1 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = predALL{newG(i,1),1}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers1{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers1{i} = tgm(1:lenG,:,:);
            end
        end
        gMarkers2 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = predALL{newG(i,1),2}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers2{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers2{i} = tgm(1:lenG,:,:);
            end
        end
        
        % transform gMarkers2 to centered coord
        for i = 1:size(gMarkers2,2)
            gm1 = gMarkers1{i};
            gm2 = gMarkers2{i};
            cm1 = gMarkers1{i}(:,1:2,5);
            cm2 = gMarkers2{i}(:,1:2,5);
            cboth = [cm1; cm2];
            mpt = median(cboth);
            
            gMarkers1{i}(:,1,:) = gMarkers1{i}(:,1,:)-mpt(1);
            gMarkers1{i}(:,2,:) = gMarkers1{i}(:,2,:)-mpt(2);
            gMarkers2{i}(:,1,:) = gMarkers2{i}(:,1,:)-mpt(1);
            gMarkers2{i}(:,2,:) = gMarkers2{i}(:,2,:)-mpt(2);
        end
        
        unitcolor = cell(4,4); newGs = reshape(isrr(newG(:,1)),[4 4]);
        for i = 1:4
            for j = 1:4
                if newGs(i,j) == 1
                    unitcolor{i,j} = [.1 .5 .1];
                else
                    unitcolor{i,j} = [.1 .5 .1];
                end

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

        allM3 = []; allJ3 = []; allC3 = []; allmc3 = [];
        for i = 1:size(newG,1)
            cx = coffx(i); cy = coffy(i); cz = coffz(i);
            NM = cat(3,gMarkers1{i}(:,:,5),gMarkers2{i}(:,:,5));
            NM(:,1,:) = NM(:,1,:)+cx;
            NM(:,2,:) = NM(:,2,:)+cy;
            NM(:,3,:) = NM(:,3,:)+cz;
            allM3 = cat(3,allM3,NM);
            newJ = [1 2]+(i-1)*2;
            allJ3 = cat(1,allJ3,newJ);
            allC3 = cat(1,allC3,[0 0 0; 0 0 0]);
        end
        skConnect.mcolor = .8*ones(32,3); skConnect.joints_idx = allJ3; skConnect.color = .8*ones(16,3);
        
        
                  close all;
        findic = figure('Name','Rat Test');
        h = cell(3,1);
        h{1} = Keypoint3DAnimator(allM3,skConnect,'Position',[0 0 1 1],'lineWidth',1.5);
        h{2} = Keypoint3DAnimator(allM,sk2,'Axes',h{1}.Axes,'lineWidth',1.5);
        h{3} = Keypoint3DAnimator(allM2,sk3,'Axes',h{1}.Axes,'lineWidth',1.5);
        
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition',[0.3934 0.0679 0.5402 0.8232]);
                view(h{1}.getAxes,20,30); axis equal
                

        x1 = 0; y1 = 0;
        for i = 1:4
            x = x1+coffx(i,1);
            for j = 1:4
                y = y1+coffy(1,j);
                r  = 350;
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + x;
                yunit = r * sin(th) + y;
                h2 = plot(xunit, yunit,'Color',unitcolor{i,j}); hold on;
            end
        end
         set(gcf,'Units','Normalized','OuterPosition',[0.3934 0.0679 0.5402 0.8232]);
                view(h{1}.getAxes,20,30); axis equal
        % axis([-200 2200 -200 2200 -200 300]);
        axis([-400 2800 -400 2800 -200 300]);
        set(gcf,'Color','white');

        savePath = ['/home/ugne/Dropbox/A_SDANNCE/bSOC_small2/' num2str(tm) '_s.avi'];
        frames = 1:lenG;

        % Uncomment to write the Animation to video.
        h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
        pause(2)
        end
    catch
    end
end






tm = 73;


G = groups40{tm};
lgs = size(G,1);
if lgs>16
    newG = G(randperm(lgs,16),:);
else
    newG = G;
end
lenG = 150;
gMarkers1 = cell(1,size(newG,1));
for i = 1:size(newG,1)
    repG = predALL{newG(i,1),1}(newG(i,2):newG(i,3),:,:);
    if size(repG,1) > lenG
        gMarkers1{i} = repG(1:lenG,:,:);
    else
        tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
        gMarkers1{i} = tgm(1:lenG,:,:);
    end
end
gMarkers2 = cell(1,size(newG,1));
for i = 1:size(newG,1)
    repG = predALL{newG(i,1),2}(newG(i,2):newG(i,3),:,:);
    if size(repG,1) > lenG
        gMarkers2{i} = repG(1:lenG,:,:);
    else
        tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
        gMarkers2{i} = tgm(1:lenG,:,:);
    end
end

mf = movieF{newG(i,1)};
vf = VideoReader(mf);
v = read(vf,[newG(i,2) newG(i,3)]);































