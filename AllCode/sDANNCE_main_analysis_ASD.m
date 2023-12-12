% ARID_L, ARID_S, ARID_SOC
% CHD8_L, CHD8_S, CHD8_SOC
% GRINB_L, GRINB_S, GRINB_SOC
% NRXN_L, NRXN_S, NRXN_SOC
% SCN2A_L, SCN2A_S, SCN2A_SOC

LONEd = [ARID_L{1}'; CHD8_L{1}'; GRINB_L{1}'; NRXN_L{1}'; SCN2A_L{1}'];
SOCd = [ARID_S{1}'; CHD8_S{1}'; GRINB_S{1}'; NRXN_S{1}'; SCN2A_S{1}'];
SSOCd = [ARID_SOC{1}'; CHD8_SOC{1}'; GRINB_SOC{1}'; NRXN_SOC{1}'; SCN2A_SOC{1}'];

LONEf = [ARID_L{2}'; CHD8_L{2}'; GRINB_L{2}'; NRXN_L{2}'; SCN2A_L{2}'];
SOCf = [ARID_S{2}'; CHD8_S{2}'; GRINB_S{2}'; NRXN_S{2}'; SCN2A_S{2}'];
SSOCf = [ARID_SOC{2}'; CHD8_SOC{2}'; GRINB_SOC{2}'; NRXN_SOC{2}'; SCN2A_SOC{2}'];

l_sig = [ARID_L(4) CHD8_L(4) GRINB_L(4) NRXN_L(4) SCN2A_L(4)];
soc_sig = [ARID_S(4) CHD8_S(4) GRINB_S(4) NRXN_S(4) SCN2A_S(4)];
ssoc_sig = [ARID_SOC(4) CHD8_SOC(4) GRINB_SOC(4) NRXN_SOC(4) SCN2A_SOC(4)];








figure(1);
imagesc(LONEd)
colormap(cmapdiff)
caxis([-.02 .02])

figure(2);
imagesc(LONEf)
colormap(cmapdiff)
caxis([-2 2])


figure(1);
imagesc(SOCd)
colormap(cmapdiff)
caxis([-.02 .02])

figure(2);
imagesc(SOCf)
colormap(cmapdiff)
caxis([-2 2])



figure(1);
imagesc(SSOCd)
colormap(cmapdiff)
caxis([-.02 .02])


figure(1);
imagesc(LONEf); set(gcf,'Position',[7 510 1249 96]); axis off
colormap(cmapdiff)
caxis([-2.5 2.5])
for group = 1:5
    for b = 1:162
        if ismember(b,l_sig{group})
             hold on; scatter(b,group,15,'k','filled')
             hold on; scatter(b,group,15,'w')
        end
    end
end



figure(2);
imagesc(SOCf); set(gcf,'Position',[7 510 1249 96]); axis off
colormap(cmapdiff)
caxis([-2.5 2.5])
for group = 1:5
    for b = 1:162
        if ismember(b,soc_sig{group})
            hold on; scatter(b,group,15,'k','filled')
             hold on; scatter(b,group,15,'w')
        end
    end
end


figure(3);
imagesc(SSOCf); set(gcf,'Position',[7 510 1249 96]); axis off
colormap(cmapdiff)
caxis([-2.5 2.5])
for group = 1:5
    for b = 1:156
        if ismember(b,ssoc_sig{group})
             hold on; scatter(b,group,15,'k','filled')
             hold on; scatter(b,group,15,'w')
        end
    end
end



d7 = [{ASD_DIFF{1,3}(:,7)} {ASD_DIFF{4,3}(:,7)} {ASD_DIFF{5,3}(:,7)}];

for b = 1:156
    figure(b)
    for i = 1:5
        y = ASD_FOLD{usestrain(i),3}(:,b);
        x = (2*i)+ones(1,length(y));
        s{i} = swarmchart(x,y,50,[0 0 0],'filled'); hold on
    end
end

% only arid, nrxn1, scn2a



b = 13;
usestrain = [1 4 5];
for i = 1:3
    y = ASD_FOLD{usestrain(i),3}(:,b); 
    x = (2*i)+ones(1,length(y));
    s{i} = swarmchart(x,y,50,[0 0 0],'filled'); hold on
end


b = 13;
usestrain = [1 4 5];
for i = 1:3
    y = ASD_DIFF{usestrain(i),3}(:,b); 
    x = (2*i)+ones(1,length(y));
    s{i} = swarmchart(x,y,50,[0 0 0],'filled'); hold on
end







b = 13;
usestrain = [1 4 5];
for i = 1:3
    y = ASD_USAGE{usestrain(i),3}(:,b);
    y2 = WT_USAGE{usestrain(i),3}(:,b);
    x = (3*i)+ones(1,length(y));
    x2 = (3*i)-1+ones(1,length(y2));
    s{i} = swarmchart(x,y,50,[0 0 0],'filled'); hold on;
    s2{i} = swarmchart(x2,y2,50,[0 0 0],'filled');
end



% social only

ASD_DIFF = cell(5,3); ASD_FOLD = cell(5,3); ASD_USAGE = cell(5,3); WT_USAGE = cell(5,3);
ASD_LFOLD = cell(5,3); WT_LFOLD = cell(5,3); WT_FOLD = cell(5,3);
for strain = 1:5
    
    x0l = find(ratGROUP==strain & ratGEN==0 & isSOC==0);
    x1l = find(ratGROUP==strain & ratGEN==1 & isSOC==0);
    
    x0s = find(ratGROUP==strain & isSOC==1 & ratGEN==0 & ratPGEN==0);
    xwt1s = find(ratGROUP==strain & isSOC==1 & ratGEN==0 & ratPGEN==1);
    xko1s = find(ratGROUP==strain & isSOC==1 & ratGEN==1 & ratPGEN==0);
    x2s = find(ratGROUP==strain & isSOC==1 & ratGEN==1 & ratPGEN==1);
    
    allwtL = histIndivF(x0l,behOrderZ)+1e-6;
    allkoL = histIndivF(x1l,behOrderZ)+1e-6;
    
    allwtS = histIndivF(x0s,behOrderZ)+1e-6;
    allwt1S = histIndivF(xwt1s,behOrderZ)+1e-6;
    allko1S = histIndivF(xko1s,behOrderZ)+1e-6;
    allkoS = histIndivF(x2s,behOrderZ)+1e-6;
    
    allwtSOC = histJointF(x0s,socOrderZ)+1e-6;
    allwt1SOC = histJointF(xwt1s,socOrderZ)+1e-6;
    allko1SOC = histJointF(xko1s,socOrderZ)+1e-6;
    allkoSOC = histJointF(x2s,socOrderZ)+1e-6;
    
    
    if strain == 1
        allwtL(13,:) = []; % remove missing/nan row
        allwt1S(24,:) = []; % remove missing/nan row
        allko1S(4,:) = [];
        allwt1SOC(24,:) = []; % remove missing/nan row
        allko1SOC(4,:) = [];
    end
    
    meanwtL = mean(allwtL);
    meanwtS = mean(allwtS);
    meanwtSOC = mean(allwtSOC);

    Ldiff = allkoL-meanwtL;
    Sdiff = allkoS-meanwtS;
    SOCdiff = allkoSOC-meanwtSOC;
    
    ASD_DIFF{strain,1} = Ldiff;
    ASD_DIFF{strain,2} = Sdiff;
    ASD_DIFF{strain,3} = SOCdiff;
    
    ASD_FOLD{strain,1} = (allkoL-meanwtL)./meanwtL;
    ASD_FOLD{strain,2} = (allkoS-meanwtS)./meanwtS;
    ASD_FOLD{strain,3} = (allkoSOC-meanwtSOC)./meanwtSOC;
    
    WT_FOLD{strain,1} = (allwtL-meanwtL)./meanwtL;
    WT_FOLD{strain,2} = (allwtS-meanwtS)./meanwtS;
    WT_FOLD{strain,3} = (allwtSOC-meanwtSOC)./meanwtSOC;

    ASD_USAGE{strain,3} = allkoSOC;
    WT_USAGE{strain,3} = allwtSOC;

    ASD_LFOLD{strain,1} = log2(allkoL./meanwtL); 
    ASD_LFOLD{strain,2} = log2(allkoS./meanwtS);
    ASD_LFOLD{strain,3} = log2(allkoSOC./meanwtSOC);

    WT_LFOLD{strain,1} = log2(allwtL./meanwtL);
    WT_LFOLD{strain,2} = log2(allkoS./meanwtS);
    WT_LFOLD{strain,3} = log2(allkoSOC./meanwtSOC);
end



allL = []; allLtag = [];
for strain = 1:5
    allL = [allL; ASD_DIFF{strain,1}]; 
    allLtag = [allLtag; strain*ones(size(ASD_DIFF{strain,1},1),1)];
end

cmapStrains = [51 102 154;
 248 181 31;
141 75 156;
11 153 73;
235 44 38]./255;


 [COEFF, SCORE, LATENT, TSQUARED,EXPLAINED] = pca(allL);

xcolors = zeros(size(allL,1),3);
for i = 1:length(allLtag)
    xcolors(i,:) = cmapStrains(allLtag(i),:);
end

h = scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),75,xcolors,'filled');
axis equal
alpha = 1;
set(h, 'MarkerFaceAlpha', alpha)
set(h, 'MarkerEdgeColor','k')



allS_SOC = []; allStag = [];
for strain = 1:5
    scombo = [ASD_DIFF{strain,2} ASD_DIFF{strain,3}];
    allS_SOC = [allS_SOC; scombo];
    allStag = [allStag; strain*ones(size(ASD_DIFF{strain,2},1),1)];
end

xcolorsS = zeros(size(allS_SOC,1),3);
for i = 1:length(allStag)
    xcolorsS(i,:) = cmapStrains(allStag(i),:);
end


 [COEFFs, SCOREs, LATENTs, TSQUAREDs,EXPLAINEDs] = pca(allS_SOC);

h = scatter3(SCOREs(:,1),SCOREs(:,2),SCOREs(:,3),100,xcolorsS,'filled');

alpha = 1;
set(h, 'MarkerFaceAlpha', alpha)
set(h, 'MarkerEdgeColor','k')

huse = find(allStag==1 | allStag==4 | allStag==5)

figure(2)
h = scatter3(SCOREs(huse,1),SCOREs(huse,2),SCOREs(huse,3),100,xcolorsS(huse,:),'filled');
axis equal
alpha = 1;
set(h, 'MarkerFaceAlpha', alpha)
set(h, 'MarkerEdgeColor','k')
axis([-0.1323    0.1062   -0.0615    0.1656   -0.1246    0.1004])
hold on


%% 


% per trial experimental changes - wtwt-koko



%% ASD MI
ASD_FC = cell(5,2); ASD_EMB = cell(5,2); 
for strain = 1:5
    
    x0l = find(ratGROUP==strain & ratGEN==0 & isSOC==0);
    x1l = find(ratGROUP==strain & ratGEN==1 & isSOC==0);
    
    x0s = find(ratGROUP==strain & isSOC==1 & ratGEN==0 & ratPGEN==0);
    xwt1s = find(ratGROUP==strain & isSOC==1 & ratGEN==0 & ratPGEN==1);
    xko1s = find(ratGROUP==strain & isSOC==1 & ratGEN==1 & ratPGEN==0);
    x2s = find(ratGROUP==strain & isSOC==1 & ratGEN==1 & ratPGEN==1);
    
    EMB_WT = cell(length(x0s),2);
    EMB_KO = cell(length(x2s),2);
    CW_WT = zeros(length(x0s),90000);
    CW_KO = zeros(length(x2s),90000);
    
    for i = 1:length(x0s)
        rec = x0s(i);
        if rec+396<1003
            rec2 = rec+396;
        else
            rec2 = rec-396;
        end
        EMB_WT{i,1} = wCOARSE{rec};
        EMB_WT{i,2} = wCOARSE{rec2};
        bbw = EMB_WT{i,1}==EMB_WT{i,2}; bb = zeros(1,90000); bb(bbw) = EMB_WT{i,1}(bbw);
        CW_WT(i,:) = bb;
    end
    
    for i = 1:length(x2s)
        rec = x2s(i);
        if rec+396<1003
            rec2 = rec+396;
        else
            rec2 = rec-396;
        end
        EMB_KO{i,1} = wCOARSE{rec};
        EMB_KO{i,2} = wCOARSE{rec2};
        bbw = EMB_KO{i,1}==EMB_KO{i,2}; bb = zeros(1,90000); bb(bbw) = EMB_KO{i,1}(bbw);
        CW_KO(i,:) = bb;
    end
    
    
    
    emb1 = combineCells(EMB_WT(:,1));
    emb2 = combineCells(EMB_WT(:,2));
    simB = [emb1 emb2];
    M = 9;
    cfmutZ = zeros(M,M);
    for i = 1:length(simB)
        currE = simB(i,:);
        idx1 = currE(1); idx2 = currE(2);
        cfmutZ(idx1,idx2) = cfmutZ(idx1,idx2)+1;
    end
    cfmutZ = cfmutZ/sum(cfmutZ(:));
    cfBG = sum(cfmutZ,2);
    
    FC_WT = zeros(M,M);
    for i = 1:M
        for j = 1:M
            FC_WT(i,j) = cfmutZ(i,j)*log2(cfmutZ(i,j)/(cfBG(i)*cfBG(j)));
        end
    end

    emb1 = combineCells(EMB_KO(:,1));
    emb2 = combineCells(EMB_KO(:,2));
    simB = [emb1 emb2];
    M = 9;
    cfmutZ = zeros(M,M);
    for i = 1:length(simB)
        currE = simB(i,:);
        idx1 = currE(1); idx2 = currE(2);
        cfmutZ(idx1,idx2) = cfmutZ(idx1,idx2)+1;
    end
    cfmutZ = cfmutZ/sum(cfmutZ(:));
    cfBG = sum(cfmutZ,2);
    
    FC_KO = zeros(M,M);
    for i = 1:M
        for j = 1:M
            FC_KO(i,j) = cfmutZ(i,j)*log2(cfmutZ(i,j)/(cfBG(i)*cfBG(j)));
        end
    end

    figure(strain); subplot(1,2,1); imagesc(FC_WT); colormap(cmapdiff); caxis([-.07 .07]);
    subplot(1,2,2); imagesc(FC_KO); colormap(cmapdiff); caxis([-.07 .07]);
    
    ASD_FC{strain,1} = FC_WT;
    ASD_FC{strain,2} = FC_KO;
    ASD_EMB{strain,1} = CW_WT;
    ASD_EMB{strain,2} = CW_KO;
    
end


for strain = 1:5
    fcwt = ASD_FC{strain,1};
    fcko = ASD_FC{strain,2};
    fcwt(isnan(fcwt)) = 0; fcko(isnan(fcwt)) = 0; fcdiff = fcko-fcwt;
    
    figure(strain); 
    subplot(1,3,1); imagesc(fcwt); colormap(cmapdiff); caxis([-.035 .035]); axis equal off
    subplot(1,3,2); imagesc(fcko); colormap(cmapdiff); caxis([-.035 .035]); axis equal off
    subplot(1,3,3); imagesc(fcdiff); colormap(cmapdiff); caxis([-.035 .035]); axis equal off
end




EMB_WT = cell(length(x0s),2);
CW_WT = zeros(length(x0s),90000);

for i = 1:length(x0s)
    rec = x0s(i);
    emb1 = ratCZB{rec,1}; emb1(emb1==0) = nan; emb1 = fillmissing(emb1,'nearest');
    emb2 = ratCZB{rec,2}; emb2(emb2==0) = nan; emb2 = fillmissing(emb2,'nearest');
    EMB_WT{i,1} = labelsCZCoarse(emb1);
    EMB_WT{i,2} = labelsCZCoarse(emb2);
    bbw = EMB_WT{i,1}==EMB_WT{i,2}; bb = zeros(1,90000); bb(bbw) = EMB_WT{i,1}(bbw);
    CW_WT(i,:) = bb; 
end

emb1 = combineCells(EMB_WT(:,1));
emb2 = combineCells(EMB_WT(:,2));
simB = [emb1 emb2];
M = 9;
cfmutZ = zeros(M,M);
for i = 1:length(simB)
    currE = simB(i,:);
    idx1 = currE(1); idx2 = currE(2);
    cfmutZ(idx1,idx2) = cfmutZ(idx1,idx2)+1;
end
cfmutZ = cfmutZ/sum(cfmutZ(:));
cfBG = sum(cfmutZ,2);

FC_WT = zeros(M,M);
for i = 1:M
    for j = 1:M
        FC_WT(i,j) = cfmutZ(i,j)*log2(cfmutZ(i,j)/(cfBG(i)*cfBG(j)));
    end
end


EMB_AMPH = cell(length(x1s),2);
CW_AMPH = zeros(length(x1s),90000);
for i = 1:length(x1s)
    rec = x1s(i);
    emb1 = ratCZB{rec,1}; emb1(emb1==0) = nan; emb1 = fillmissing(emb1,'nearest');
    emb2 = ratCZB{rec,2}; emb2(emb2==0) = nan; emb2 = fillmissing(emb2,'nearest');
    EMB_AMPH{i,1} = labelsCZCoarse(emb1);
    EMB_AMPH{i,2} = labelsCZCoarse(emb2);
    bbw = EMB_AMPH{i,1}==EMB_AMPH{i,2}; bb = zeros(1,90000); bb(bbw) = EMB_AMPH{i,1}(bbw);
    CW_AMPH(i,:) = bb; 
end

emb1 = combineCells(EMB_AMPH(:,1));
emb2 = combineCells(EMB_AMPH(:,2));
simB = [emb1 emb2];
M = 9;
cfmutZ = zeros(M,M);
for i = 1:length(simB)
    currE = simB(i,:);
    idx1 = currE(1); idx2 = currE(2);
    cfmutZ(idx1,idx2) = cfmutZ(idx1,idx2)+1;
end
cfmutZ = cfmutZ/sum(cfmutZ(:));
cfBG2 = sum(cfmutZ,2);
cfBG1 = sum(cfmutZ,1);

FC_AMPH = zeros(M,M);
for i = 1:M
    for j = 1:M
        FC_AMPH(i,j) = cfmutZ(i,j)*log2(cfmutZ(i,j)/(cfBG2(i)*cfBG1(j)));
    end
end


EMB_AMPH2 = cell(length(x2s),2);
CW_AMPH2 = zeros(length(x2),90000);
for i = 1:length(x2s)
    rec = x2s(i);
    emb1 = ratCZB{rec,1}; emb1(emb1==0) = nan; emb1 = fillmissing(emb1,'nearest');
    emb2 = ratCZB{rec,2}; emb2(emb2==0) = nan; emb2 = fillmissing(emb2,'nearest');
    EMB_AMPH2{i,1} = labelsCZCoarse(emb1);
    EMB_AMPH2{i,2} = labelsCZCoarse(emb2);
    bbw = EMB_AMPH2{i,1}==EMB_AMPH2{i,2}; bb = zeros(1,90000); bb(bbw) = EMB_AMPH2{i,1}(bbw);
    CW_AMPH2(i,:) = bb;
end

emb1 = combineCells(EMB_AMPH2(:,1));
emb2 = combineCells(EMB_AMPH2(:,2));
simB = [emb1 emb2];
M = 9;
cfmutZ = zeros(M,M);
for i = 1:length(simB)
    currE = simB(i,:);
    idx1 = currE(1); idx2 = currE(2);
    cfmutZ(idx1,idx2) = cfmutZ(idx1,idx2)+1;
end
cfmutZ = cfmutZ/sum(cfmutZ(:));
cfBG2 = sum(cfmutZ,2);
cfBG1 = sum(cfmutZ,1);

FC_AMPH2 = zeros(M,M);
for i = 1:M
    for j = 1:M
        FC_AMPH2(i,j) = cfmutZ(i,j)*log2(cfmutZ(i,j)/(cfBG2(i)*cfBG1(j)));
    end
end

figure(1); imagesc(FC_WT); colormap(cmapdiff); caxis([-.07 .07]); axis equal off
saveas(gcf,'A_SDANNCE/bigRun/figs0328/MI_WT.tif'); 
colorbar
saveas(gcf,'A_SDANNCE/bigRun/figs0328/MI_WT_CB.tif'); 

figure(2); imagesc(FC_AMPH); colormap(cmapdiff); caxis([-.07 .07]); axis equal off
saveas(gcf,'A_SDANNCE/bigRun/figs0328/MI_AMPH.tif'); 
colorbar
saveas(gcf,'A_SDANNCE/bigRun/figs0328/MI_AMPH_CB.tif'); 






















