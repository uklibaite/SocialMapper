%% some basic metrics
%212mm radius for approximately half the area
% distance traveled, time spent in center
DT = nan(1002,1); DC = nan(1002,1);
for i = 1:1002
    try
    pred = predALL{i,1};
    c1 = squeeze(pred(:,:,5)); c1s = medfilt1(c1,50);
    xmin = prctile(c1(:,1),1); xmax = prctile(c1(:,1),99);
    ymin = prctile(c1(:,2),1); ymax = prctile(c1(:,2),99);
    xcent = (xmax+xmin)./2; ycent = (ymax+ymin)./2;
    cdiff = [c1(:,1)-xcent c1(:,2)-ycent];
    dc = sqrt(sum(cdiff.^2,2));
    dcfrac = mean(dc<212);
    trav = sum(sqrt(diff(c1s(:,[1 2])).^2),2);
    dt = sum(trav(2:end));
    DT(i) = dt;
    DC(i) = dcfrac;
    catch
    end
end

% fraction of time spent in 3 proximity bins
PROXFRAC = nan(1002,3);
for i = 1:396
    dc = allDC{i};
    f1 = sum(dc<200)./length(dc); 
    f2 = sum(dc>200 & dc<400)./length(dc);
    f3 = sum(dc>400)./length(dc);
    PROXFRAC(210+i,:) = [f1 f2 f3];
    PROXFRAC(210+i+396,:) = [f1 f2 f3];
end


% find wt/amphetamine groups
x0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0);
x0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1);
x1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0);
x1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1);
x2l = find(ratGROUP==6 & isAMPH==2 & isSOC==0);
x2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1);


% distance traveled for 5 wt/amphetamine groups - swarm plot
dt0l = DT(x0l);
dt1l = DT(x1l);
dt0s = DT(x0s);
dt1s = DT(x1s);
dt2s = DT(x2s);

x1 = ones(1,length(x0l));
x2 = 2+ones(1,length(x1l));
x3 = 4+ones(1,length(x0s));
x4 = 6+ones(1,length(x1s));
x5 = 8+ones(1,length(x2s));
s1 = swarmchart(x1,dt0l,50,[.3 .3 .3],'filled'); hold on
s2 = swarmchart(x2,dt1l,50,[0 1 1],'filled');
s3 = swarmchart(x3,dt0s,50,[0 0 1],'filled');
s4 = swarmchart(x4,dt1s,50,[1 0 0],'filled');
s5 = swarmchart(x5,dt2s,50,[0 1 0],'filled');

s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];
s3.MarkerEdgeColor = [0 0 0];
s4.MarkerEdgeColor = [0 0 0];
s5.MarkerEdgeColor = [0 0 0];

% 
figure(2)
dt0l = DC(x0l);
dt1l = DC(x1l);
dt0s = DC(x0s);
dt1s = DC(x1s);
dt2s = DC(x2s);

x1 = ones(1,length(x0l));
x2 = 2+ones(1,length(x1l));
x3 = 4+ones(1,length(x0s));
x4 = 6+ones(1,length(x1s));
x5 = 8+ones(1,length(x2s));
s1 = swarmchart(x1,dt0l,50,[.3 .3 .3],'filled'); hold on
s2 = swarmchart(x2,dt1l,50,[0 1 1],'filled');
s3 = swarmchart(x3,dt0s,50,[0 0 1],'filled');
s4 = swarmchart(x4,dt1s,50,[1 0 0],'filled');
s5 = swarmchart(x5,dt2s,50,[0 1 0],'filled');

s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];
s3.MarkerEdgeColor = [0 0 0];
s4.MarkerEdgeColor = [0 0 0];
s5.MarkerEdgeColor = [0 0 0];

y1 = allwtL(:,100);
y2 = allamphL(:,100);

x1 = ones(1,24);
x2 = 2+ones(1,6);

s1 = swarmchart(x1,y1,50,[.5 .5 .5],'filled'); hold on
s2 = swarmchart(x2,y2,50,[.5 .5 .5],'filled');
s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];


x1 = ones(1,24); y1 = PROXFRAC(x0s(1:24),1);
x2 = 2+ones(1,21); y2 = PROXFRAC(x1s(:),1);
x3 = 5+ones(1,24); y3 = PROXFRAC(x0s(1:24),2);
x4 = 7+ones(1,21); y4 = PROXFRAC(x1s(:),2);
x5 = 10+ones(1,24); y5 = PROXFRAC(x0s(1:24),3);
x6 = 12+ones(1,21); y6 = PROXFRAC(x1s(:),3);
s1 = swarmchart(x1,y1,50,[.3 .3 .3],'filled'); hold on;
s2 = swarmchart(x2,y2,50,[.3 .3 .3],'filled'); hold on;
s3 = swarmchart(x3,y3,50,[.3 .3 .3],'filled'); hold on;
s4 = swarmchart(x4,y4,50,[.3 .3 .3],'filled'); hold on;
s5 = swarmchart(x5,y5,50,[.3 .3 .3],'filled'); hold on;
s6 = swarmchart(x6,y6,50,[.3 .3 .3],'filled'); hold on;
s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];
s3.MarkerEdgeColor = [0 0 0];
s4.MarkerEdgeColor = [0 0 0];
s5.MarkerEdgeColor = [0 0 0];
s6.MarkerEdgeColor = [0 0 0];



wtIndivL = histIndivC(x0l,:); MwtIndivL = mean(wtIndivL);
amphIndivL = histIndivC(x1l,:); MamphIndivL = mean(amphIndivL);
wtIndivS = histIndivC(x0s,:); MwtIndivS = mean(wtIndivS);
amphIndivS = histIndivC(x1s,:); MamphIndivS = mean(amphIndivS);
partIndivS = histIndivC(x2s,:); MpartIndivS = mean(partIndivS);

wtJointSOC = histJointC(x0s,:); MwtJointSOC = mean(wtJointSOC);
amphJointSOC = histJointC(x1s,:); MamphJointSOC = mean(amphJointSOC);
partJointSOC = histJointC(x2s,:); MpartJointSOC = mean(partJointSOC);





% diff1 = WT-Lone compared to AMPH-Lone and WT-Social
% diff2 = WT-Social compared to AMPH-Social and PARTNER-Social
diff1 = zeros(9,2);
for b = 1:9
    lb = MwtIndivL(b);
    lbamph = MamphIndivL(b);
    lbsoc = MwtIndivS(b);
    damph = (lbamph-lb)./lb;
    dsoc = (lbsoc-lb)./lb;
    diff1(b,:) = [damph dsoc];
end


diff2 = zeros(9,2);
for b = 1:9
    lb = MwtIndivS(b);
    lbamph = MamphIndivS(b);
    lbsoc = MpartIndivS(b);
    damph = (lbamph-lb)./lb;
    dsoc = (lbsoc-lb)./lb;
    diff2(b,:) = [damph dsoc];
end


diff3 = zeros(9,2);
for b = 1:9
    lb = MwtJointSOC(b);
    lbamph = MamphJointSOC(b);
    lbpart = MpartJointSOC(b);
    damph = (lbamph-lb)./lb;
    dpart = (lbpart-lb)./lb;
    diff3(b,:) = [damph dpart];
end






r = 1:27;
dr = reshape(r,[3,9])';
dr2 = dr(:,1:2);


b = bar(dr2(:), diff1(:))
b.FaceColor = 'flat';
b.CData(1:2:end,:) = repmat([0 1 1],[9 1]);
b.CData(2:2:end,:) = repmat([0 0 1],[9 1]);

axis([-.2 27.2 -2 2])
set(gcf,'Position',[340 331 1261 463]);
saveas(gcf,'figs20230917/LE_diff1.tif');
axis off
saveas(gcf,'figs20230917/LE_diff1_NA.tif');


b = bar(dr2(:), diff2(:))
b.FaceColor = 'flat';
b.CData(1:2:end,:) = repmat([1 0 0],[9 1]);
b.CData(2:2:end,:) = repmat([0 1 0],[9 1]);

axis([-.2 27.2 -1 1])
set(gcf,'Position',[340 331 1261 463]);
saveas(gcf,'figs20230917/LE_SOC_diff2.tif');
axis off
saveas(gcf,'figs20230917/LE_SOC_diff2_NA.tif');



r2 = 1:63;
dr = reshape(r2,[7,9])';
dr3 = dr(:,[1 3]);


bx = dr3(1:7,:);
by = diff3(btags,:);
b = bar(bx(:), by(:));
b.FaceColor = 'flat';
b.CData(1:2:end,:) = repmat([1 1 1],[7 1]);
b.CData(2:2:end,:) = repmat([1 1 1],[7 1]);

axis([-.2 46.2 -3 3])
set(gcf,'Position',[340 331 1261 463]);
saveas(gcf,'figs20230917/LE_SOC_diff3BW.tif');
axis off
saveas(gcf,'figs20230917/LE_SOC_diff3BW_NA.tif');



%% BAR+SWARM
btags = [1 4 2 5 3 8 9];
hold on

allAMPH = []; allPART = [];
for rn = 1:6
    rnx = rn-1;
    xi0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0 & ratIDN==rnx);
    xi1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0 & ratIDN==rnx);
    xi0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1 & ratIDN==rnx);
    xi1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1 & ratIDN==rnx);
    xi2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1 & ratIDN==rnx);
    
    wtJC = histJointC(xi0s,btags); mwtJC = mean(wtJC);
    amphJC = histJointC(xi1s,btags); partJC = histJointC(xi2s,btags);
    damphJC = amphJC-mwtJC; dpartJC = partJC-mwtJC;
    %damphJC = log2(amphJC./mwtJC); dpartJC = log2(partJC./mwtJC);
    allAMPH = [allAMPH; damphJC];
    allPART = [allPART; dpartJC];
end

C = cell(2,7);
for i = 1:7
    C{1,i} = allAMPH(:,i);
    C{2,i} = allPART(:,i);
end

cents = [1 2 4 5 7 8 10 11 13 14 16 17 19 20];
cc = repmat([0 0 0; .5 .5 .5],[7 1]);
[outinfo,tt] = mySwarmNA(C(:),cents,.5,cc);
axis([0 21 -.4 .4])
axis([0 21 -3.5 3.5])


for i = 1:7
    figure(i)
    cents = [1 2];
    cc = [0 0 0; .5 .5 .5];
    [outinfo,tt] = mySwarmNA(C(:,i),cents,.5,cc);
    max1 = max(combineCells(C(:,i)))*1.1;
    min1 = min(combineCells(C(:,i)))*1.1;
    axis([0 3 min1 max1]);
    set(gcf,'Position',[1377 83 119 420])
end



for ib = 1:7
    b = btags(ib);
    wtb = wtJointSOC(:,b); mwtb = MwtJointSOC(b);
    amphb = amphJointSOC(:,b);
    partb = partJointSOC(:,b);
    
    xb = ones(1,length(amphb))+dr3(ib,1)-1;
    yb = (amphb-mwtb)./mwtb;
    
    xb2 = ones(1,length(partb))+dr3(ib,2)-1;
    yb2 = (partb-mwtb)./mwtb;
    
    s{b,1} = swarmchart(xb,yb,10,[0 0 0],'filled'); 
    s{b,2} = swarmchart(xb2,yb2,10,[0 0 0],'filled'); hold on
end


%% sDANNCE_main analysis 
histIndivF = nan(1002,163); histIndivC = nan(1002,9);
histJointF = nan(1002,161); histJointC = nan(1002,9);
histSOCC = nan(1002,4);

for i = 1:1002
    
    indivF = wFINE{i};
    indivC = wCOARSE{i};
    jointF = wsocFINE{i};
    jointC = wsocCOARSE{i};
    socC = SMSOC{i};
    hif = hist(indivF,1:163);
    hic = hist(indivC,1:9);
    hjf = hist(jointF,1:161);
    hjc = hist(jointC,0:8);
    hsc = hist(socC,1:4);
    try
        histIndivF(i,:) = hif./sum(hif);
        histIndivC(i,:) = hic./sum(hic);
    catch
    end
    try
        histJointF(i,:) = hjf./sum(hjf);
        histJointC(i,:) = hjc./sum(hjc);
        histSOCC(i,:) = hsc./sum(hsc);
    catch
    end
end

% csv with labels 

T = readtable('/home/ugne/Downloads/individualClusterDescriptions - Sheet1.csv');
indIndex = table2array(T(1:163,2));
desc = T(2:163,3);
dnum = T(2:163,2);
indText = cell(162,1);
for i = 1:162
    indText{i} = [num2str(table2array(dnum(i,1))) '. ' desc{i,1}{1}];
end

TSOC = readtable('/home/ugne/Downloads/socialClusterDescriptions - Sheet1.csv');
jointIndex = table2array(TSOC(1:156,2));
desc = TSOC(:,3);
dnum = TSOC(:,2);
jointText = cell(156,1);
for i = 1:156
    jointText{i} = [num2str(table2array(dnum(i,1))) '. ' desc{i,1}{1}];
end

colorsInd = zeros(162,3);
colorsJoint = zeros(156,3);
for i = 1:length(colorsInd)
    colorsInd(i,:) = colorsCoarse(behOrderTagZ(i)+1,:);
end
for i = 1:length(colorsJoint)
    colorsJoint(i,:) = socColors(socOrderTagZ(i)+2,:);
end


x0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0);
x0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1);
x1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0);
x1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1);
x2l = find(ratGROUP==6 & isAMPH==2 & isSOC==0);
x2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1);

allwtL = histIndivF(x0l,behOrderZ)+1e-6;
allamphL = histIndivF(x1l,behOrderZ)+1e-6;
allpartL = histIndivF(x2l,behOrderZ)+1e-6;

allwtS = histIndivF(x0s,behOrderZ)+1e-6;
allamphS = histIndivF(x1s,behOrderZ)+1e-6;
allpartS = histIndivF(x2s,behOrderZ)+1e-6;

allwtSOC = histJointF(x0s,socOrderZ)+1e-6;
allamphSOC = histJointF(x1s,socOrderZ)+1e-6;
allpartSOC = histJointF(x2s,socOrderZ)+1e-6;

nnum = 100;
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allamphL, nnum);
WT_AMPH_L = [{mdiff} {mfold} {fdrm} {idc}];

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allamphS, nnum);
WT_AMPH_S = [{mdiff} {mfold} {fdrm} {idc}]; 

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allamphSOC, nnum);
WT_AMPH_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 

figure(1); dispBehaviorSet(WT_AMPH_L,colorsInd,indText);
axis off
figure(2); dispBehaviorSet(WT_AMPH_S,colorsInd,indText);
axis off
figure(3); dispBehaviorSet(WT_AMPH_SOC,colorsJoint,jointText);
axis off




nnum = 100;
% ARID
x0l = find(ratGROUP==1 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==1 & ratGEN==1 & isSOC==0);
x0s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ)+1e-6; allwtL(13,:) = []; % remove missing/nan row
allkoL = histIndivF(x1l,behOrderZ)+1e-6;

allwtS = histIndivF(x0s,behOrderZ)+1e-6;
allwt1S = histIndivF(xwt1s,behOrderZ)+1e-6; allwt1S(24,:) = []; % remove missing/nan row
allko1S = histIndivF(xko1s,behOrderZ)+1e-6; allko1S(4,:) = [];
allkoS = histIndivF(x2s,behOrderZ)+1e-6;

allwtSOC = histJointF(x0s,socOrderZ)+1e-6;
allwt1SOC = histJointF(xwt1s,socOrderZ)+1e-6; allwt1SOC(24,:) = []; % remove missing/nan row
allko1SOC = histJointF(xko1s,socOrderZ)+1e-6; allko1SOC(4,:) = [];
allkoSOC = histJointF(x2s,socOrderZ)+1e-6;

[mdiffL, mfoldL, fdrmL, idcL, nupL, ndownL] = testBehaviorSet5(allwtL, allkoL ,nnum);
[mdiffS, mfoldS, fdrmS, idcS, nupS, ndownS] = testBehaviorSet5(allwtS, allkoS ,nnum);
[mdiffSOC, mfoldSOC, fdrmSOC, idcSOC, nupSOC, ndownSOC] = testBehaviorSet5(allwtSOC, allkoSOC ,nnum);

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allkoL ,nnum); ARID_L = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allkoS ,nnum); ARID_S = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allkoSOC ,nnum); ARID_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 

% figure(1); dispBehaviorSet(ARID_L,colorsInd,indText); axis off; 
% figure(2); dispBehaviorSet(ARID_S,colorsInd,indText); axis off;
% figure(4); dispBehaviorSet(ARID_SOC,colorsJoint,jointText); axis off;


% CHD8
x0l = find(ratGROUP==2 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==2 & ratGEN==1 & isSOC==0);
x0s = find(ratGROUP==2 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==2 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==2 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==2 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

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

[mdiffL, mfoldL, fdrmL, idcL, nupL, ndownL] = testBehaviorSet5(allwtL, allkoL ,nnum);
[mdiffS, mfoldS, fdrmS, idcS, nupS, ndownS] = testBehaviorSet5(allwtS, allkoS ,nnum);
[mdiffSOC, mfoldSOC, fdrmSOC, idcSOC, nupSOC, ndownSOC] = testBehaviorSet5(allwtSOC, allkoSOC ,nnum);

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allkoL ,nnum); CHD8_L = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allkoS ,nnum); CHD8_S = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allkoSOC ,nnum); CHD8_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 

% figure(1); dispBehaviorSet(CHD8_L,colorsInd,indText); axis off; 
% figure(2); dispBehaviorSet(CHD8_S,colorsInd,indText); axis off; 
% figure(3); dispBehaviorSet(CHD8_SOC,colorsJoint,jointText); axis off; 

% GRINB
x0l = find(ratGROUP==3 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==3 & ratGEN==1 & isSOC==0);
x0s = find(ratGROUP==3 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==3 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==3 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==3 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

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


[mdiffL, mfoldL, fdrmL, idcL, nupL, ndownL] = testBehaviorSet5(allwtL, allkoL ,nnum);
[mdiffS, mfoldS, fdrmS, idcS, nupS, ndownS] = testBehaviorSet5(allwtS, allkoS ,nnum);
[mdiffSOC, mfoldSOC, fdrmSOC, idcSOC, nupSOC, ndownSOC] = testBehaviorSet5(allwtSOC, allkoSOC ,nnum);

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allkoL ,nnum); GRINB_L = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allkoS ,nnum); GRINB_S = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allkoSOC ,nnum); GRINB_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 

% figure(1); dispBehaviorSet(GRINB_L,colorsInd,indText); axis off; 
% figure(2); dispBehaviorSet(GRINB_S,colorsInd,indText); axis off; 
% figure(3); dispBehaviorSet(GRINB_SOC,colorsJoint,jointText); axis off; 

% NRXN1
x0l = find(ratGROUP==4 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==4 & ratGEN==1 & isSOC==0);
x0s = find(ratGROUP==4 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==4 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==4 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==4 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

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


[mdiffL, mfoldL, fdrmL, idcL, nupL, ndownL] = testBehaviorSet5(allwtL, allkoL ,nnum);
[mdiffS, mfoldS, fdrmS, idcS, nupS, ndownS] = testBehaviorSet5(allwtS, allkoS ,nnum);
[mdiffSOC, mfoldSOC, fdrmSOC, idcSOC, nupSOC, ndownSOC] = testBehaviorSet5(allwtSOC, allkoSOC ,nnum);

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allkoL ,nnum); NRXN_L = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allkoS ,nnum); NRXN_S = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allkoSOC ,nnum); NRXN_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 


% figure(1); dispBehaviorSet(NRXN_L,colorsInd,indText); axis off; 
% figure(2); dispBehaviorSet(NRXN_S,colorsInd,indText); axis off; 
% figure(3); dispBehaviorSet(NRXN_SOC,colorsJoint,jointText); axis off; 

% SCN2A
x0l = find(ratGROUP==5 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==5 & ratGEN==1 & isSOC==0);
x0s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

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


[mdiffL, mfoldL, fdrmL, idcL, nupL, ndownL] = testBehaviorSet5(allwtL, allkoL ,nnum);
[mdiffS, mfoldS, fdrmS, idcS, nupS, ndownS] = testBehaviorSet5(allwtS, allkoS ,nnum);
[mdiffSOC, mfoldSOC, fdrmSOC, idcSOC, nupSOC, ndownSOC] = testBehaviorSet5(allwtSOC, allkoSOC ,nnum);

[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtL, allkoL ,nnum); SCN2A_L = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtS, allkoS ,nnum); SCN2A_S = [{mdiff} {mfold} {fdrm} {idc}]; 
[mdiff, mfold, fdrm, idc] = testBehaviorSet(allwtSOC, allkoSOC ,nnum); SCN2A_SOC = [{mdiff} {mfold} {fdrm} {idc}]; 

% figure(1); dispBehaviorSet(SCN2A_L,colorsInd,indText); axis off; 
% figure(2); dispBehaviorSet(SCN2A_S,colorsInd,indText); axis off; 
% figure(3); dispBehaviorSet(SCN2A_SOC,colorsJoint,jointText); axis off; 






b = bar(SCN2A_L{2}); hold on;
b.FaceColor = 'flat';
for i = 1:162
    b.CData(i,:) = colorsCoarse(behOrderTagZ(i)+1,:);
end

figure(1); imagesc(mdiff'); colormap(cmapdiff); caxis([-.05 .05])
figure(2); imagesc(mfold'); colormap(cmapdiff); caxis([-3 3])




bc = bar(mdiff(idc)); hold on;
bc.FaceColor = 'flat';
for i = 1:length(idc)
    bc.CData(i,:) = colorsCoarse(behOrderTagZ(idc(i))+1,:);
end

for i = 1:length(idc)
    if mdiff(idc(i))<0
        t = text(i,-.19,indText(idc(i)));
        t.Rotation = 90; t.FontSize = 7;
    else
        t = text(i,.05,indText(idc(i)));
        t.Rotation = 90; t.FontSize = 7;
    end
end
axis([-.2 41.2 -.2 .2])
set(gcf,'Position',[364 68 1557 952]);


    %C = [{G1} {G2}]; cents = [1 3];
    %[outinfo,tt] = mySwarmNA(C,cents,.5,repmat([0 0 0],[2 1]));
    
    %figure(2)


%% WT-AMPHETAMINE-PARTNER DATA
x0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0);
x0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1);
x1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0);
x1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1);
x2l = find(ratGROUP==6 & isAMPH==2 & isSOC==0);
x2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1);

allwtL = histIndivF(x0l,behOrderZ)+1e-6;
allamphL = histIndivF(x1l,behOrderZ)+1e-6;
allpartL = histIndivF(x2l,behOrderZ)+1e-6;

allwtS = histIndivF(x0s,behOrderZ)+1e-6;
allamphS = histIndivF(x1s,behOrderZ)+1e-6;
allpartS = histIndivF(x2s,behOrderZ)+1e-6;

allwtSOC = histJointF(x0s,socOrderZ)+1e-6;
allamphSOC = histJointF(x1s,socOrderZ)+1e-6;
allpartSOC = histJointF(x2s,socOrderZ)+1e-6;


nnum = 1000000;
[pvalue, testStat, allts] = testGroup(allwtL,allamphL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.12742 numGreater = 0


[pvalue, testStat, allts] = testGroup(allwtL,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])




[pvalue, testStat, allts] = testGroup(allwtS,allamphS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.1019 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtS,allpartS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.011807 numGreater = 383

[pvalue, testStat, allts] = testGroup(allamphS,allpartS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.10301 numGreater = 0


[pvalue, testStat, allts] = testGroup(allwtSOC,allamphSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.049846 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtSOC,allpartSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.054382 numGreater = 0

[pvalue, testStat, allts] = testGroup(allamphSOC,allpartSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.093552 numGreater = 0












%% ARID1B
nnum = 1000000;

x0l = find(ratGROUP==1 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==1 & ratGEN==1 & isSOC==0);

x0s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ); allwtL(13,:) = []; % remove missing/nan row
allkoL = histIndivF(x1l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allwt1S = histIndivF(xwt1s,behOrderZ); allwt1S(24,:) = []; % remove missing/nan row
allko1S = histIndivF(xko1s,behOrderZ); allko1S(4,:) = [];
allkoS = histIndivF(x2s,behOrderZ);

allwtSOC = histJointF(x0s,socOrderZ);
allwt1SOC = histJointF(xwt1s,socOrderZ); allwt1SOC(24,:) = []; % remove missing/nan row
allko1SOC = histJointF(xko1s,socOrderZ); allko1SOC(4,:) = [];
allkoSOC = histJointF(x2s,socOrderZ);


[pvalue, testStat, allts] = testGroup(allwtL,allkoL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.024594 numGreater = 122


[pvalue, testStat, allts] = testGroup(allwtS,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.034676 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtS,allwt1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0035285 numGreater = 111195

[pvalue, testStat, allts] = testGroup(allko1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.001884 numGreater = 793413

[pvalue, testStat, allts] = testGroup(allwt1S,allko1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.024338 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.027112 numGreater = 0

[pvalue, testStat, allts] = testGroup(allko1S,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.030296 numGreater = 0




[pvalue, testStat, allts] = testGroup(allwtSOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.043481 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtSOC,allwt1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.024188 numGreater = 0

[pvalue, testStat, allts] = testGroup(allko1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.01987 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1SOC,allko1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.032597 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.019964 numGreater = 0

[pvalue, testStat, allts] = testGroup(allko1SOC,allwtSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.023466 numGreater = 0


%% CHD8
nnum = 1000000;

x0l = find(ratGROUP==2 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==2 & ratGEN==1 & isSOC==0);

x0s = find(ratGROUP==2 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==2 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==2 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==2 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ); 
allkoL = histIndivF(x1l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allwt1S = histIndivF(xwt1s,behOrderZ);
allko1S = histIndivF(xko1s,behOrderZ); 
allkoS = histIndivF(x2s,behOrderZ);

allwtSOC = histJointF(x0s,socOrderZ);
allwt1SOC = histJointF(xwt1s,socOrderZ); 
allko1SOC = histJointF(xko1s,socOrderZ); 
allkoSOC = histJointF(x2s,socOrderZ);


[pvalue, testStat, allts] = testGroup(allwtL,allkoL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0086734 numGreater = 37430


[pvalue, testStat, allts] = testGroup(allwtS,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.017372 numGreater = 2

[pvalue, testStat, allts] = testGroup(allwtS,allwt1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0027271 numGreater = 132043

[pvalue, testStat, allts] = testGroup(allko1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0013128 numGreater = 847871

[pvalue, testStat, allts] = testGroup(allwt1S,allko1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.015829 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.012556 numGreater = 2

[pvalue, testStat, allts] = testGroup(allko1S,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.020731 numGreater = 0




[pvalue, testStat, allts] = testGroup(allwtSOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.017757 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtSOC,allwt1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.010287 numGreater = 136

[pvalue, testStat, allts] = testGroup(allko1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.010191 numGreater = 806

[pvalue, testStat, allts] = testGroup(allwt1SOC,allko1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.015585 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0097058 numGreater = 1563

[pvalue, testStat, allts] = testGroup(allko1SOC,allwtSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.01084 numGreater = 50

%% GRINB
nnum = 1000000;

x0l = find(ratGROUP==3 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==3 & ratGEN==1 & isSOC==0);

x0s = find(ratGROUP==3 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==3 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==3 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==3 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ); 
allkoL = histIndivF(x1l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allwt1S = histIndivF(xwt1s,behOrderZ);
allko1S = histIndivF(xko1s,behOrderZ); 
allkoS = histIndivF(x2s,behOrderZ);

allwtSOC = histJointF(x0s,socOrderZ);
allwt1SOC = histJointF(xwt1s,socOrderZ); 
allko1SOC = histJointF(xko1s,socOrderZ); 
allkoSOC = histJointF(x2s,socOrderZ);


[pvalue, testStat, allts] = testGroup(allwtL,allkoL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.014339 numGreater = 39734


[pvalue, testStat, allts] = testGroup(allwtS,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.017615 numGreater = 12

[pvalue, testStat, allts] = testGroup(allwtS,allwt1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0018011 numGreater = 853983

[pvalue, testStat, allts] = testGroup(allko1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0017726 numGreater = 808048

[pvalue, testStat, allts] = testGroup(allwt1S,allko1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.014793 numGreater = 2

[pvalue, testStat, allts] = testGroup(allwt1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.017268 numGreater = 1

[pvalue, testStat, allts] = testGroup(allko1S,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.015686 numGreater = 16




[pvalue, testStat, allts] = testGroup(allwtSOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.01107 numGreater = 50

[pvalue, testStat, allts] = testGroup(allwtSOC,allwt1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0070782 numGreater = 33406

[pvalue, testStat, allts] = testGroup(allko1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0063993 numGreater = 43914

[pvalue, testStat, allts] = testGroup(allwt1SOC,allko1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0071899 numGreater = 6885

[pvalue, testStat, allts] = testGroup(allwt1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0070523 numGreater = 21844

[pvalue, testStat, allts] = testGroup(allko1SOC,allwtSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0073176 numGreater = 22308


%% NRXN1
nnum = 1000000;

x0l = find(ratGROUP==4 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==4 & ratGEN==1 & isSOC==0);

x0s = find(ratGROUP==4 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==4 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==4 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==4 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ); 
allkoL = histIndivF(x1l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allwt1S = histIndivF(xwt1s,behOrderZ);
allko1S = histIndivF(xko1s,behOrderZ); 
allkoS = histIndivF(x2s,behOrderZ);

allwtSOC = histJointF(x0s,socOrderZ);
allwt1SOC = histJointF(xwt1s,socOrderZ); 
allko1SOC = histJointF(xko1s,socOrderZ); 
allkoSOC = histJointF(x2s,socOrderZ);


[pvalue, testStat, allts] = testGroup(allwtL,allkoL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.028237 numGreater = 1


[pvalue, testStat, allts] = testGroup(allwtS,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.024934 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtS,allwt1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.002036 numGreater = 246177

[pvalue, testStat, allts] = testGroup(allko1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0021818 numGreater = 350060

[pvalue, testStat, allts] = testGroup(allwt1S,allko1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.013722 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.020037 numGreater = 0

[pvalue, testStat, allts] = testGroup(allko1S,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.016954 numGreater = 0




[pvalue, testStat, allts] = testGroup(allwtSOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.022996 numGreater = 1

[pvalue, testStat, allts] = testGroup(allwtSOC,allwt1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.010509 numGreater = 287

[pvalue, testStat, allts] = testGroup(allko1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0083033 numGreater = 12393

[pvalue, testStat, allts] = testGroup(allwt1SOC,allko1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0067751 numGreater = 10593

[pvalue, testStat, allts] = testGroup(allwt1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0079308 numGreater = 15370

[pvalue, testStat, allts] = testGroup(allko1SOC,allwtSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.010155 numGreater = 533

%% SCN2A
nnum = 1000000;

x0l = find(ratGROUP==5 & ratGEN==0 & isSOC==0);
x1l = find(ratGROUP==5 & ratGEN==1 & isSOC==0);

x0s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==0);
xwt1s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==1); 
xko1s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==0); 
x2s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==1); 

allwtL = histIndivF(x0l,behOrderZ); 
allkoL = histIndivF(x1l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allwt1S = histIndivF(xwt1s,behOrderZ);
allko1S = histIndivF(xko1s,behOrderZ); 
allkoS = histIndivF(x2s,behOrderZ);

allwtSOC = histJointF(x0s,socOrderZ);
allwt1SOC = histJointF(xwt1s,socOrderZ); 
allko1SOC = histJointF(xko1s,socOrderZ); 
allkoSOC = histJointF(x2s,socOrderZ);


[pvalue, testStat, allts] = testGroup(allwtL,allkoL,nnum);
plot(1:20:length(allts),allts(1:20:end),'LineWidth',.06,'Color','k')
hold on; plot([0 nnum],[testStat testStat],'--r','LineWidth',2)
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.053716 numGreater = 67


[pvalue, testStat, allts] = testGroup(allwtS,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.02141 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtS,allwt1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0043377 numGreater = 112127

[pvalue, testStat, allts] = testGroup(allko1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.0027797 numGreater = 364235

[pvalue, testStat, allts] = testGroup(allwt1S,allko1S,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.037026 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1S,allkoS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.032577 numGreater = 0

[pvalue, testStat, allts] = testGroup(allko1S,allwtS,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.024671 numGreater = 0




[pvalue, testStat, allts] = testGroup(allwtSOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.029792 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwtSOC,allwt1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.018001 numGreater = 2

[pvalue, testStat, allts] = testGroup(allko1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.016946 numGreater = 21

[pvalue, testStat, allts] = testGroup(allwt1SOC,allko1SOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.026766 numGreater = 0

[pvalue, testStat, allts] = testGroup(allwt1SOC,allkoSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.016884 numGreater = 33

[pvalue, testStat, allts] = testGroup(allko1SOC,allwtSOC,nnum);
fprintf(['testStat = ' num2str(testStat) ' numGreater = ' num2str(sum(allts>testStat)) '\n'])
% testStat = 0.017558 numGreater = 2











nsamps = 10000;
mwtL = mean(allwtL,'omitnan');
mamphL = mean(allamphL,'omitnan');
mpartL = mean(allpartL,'omitnan');

wtL_comp = zeros(162,5);
for b = 1:162
    [m1,m2,xt] = findBootMean(allwtL(:,b),allamphL(:,b),nsamps);
    wtL_comp(b,1) = m1; wtL_comp(b,2) = m2; wtL_comp(b,3) = xt./nsamps;
end
pvalues = mattest(allwtL',allamphL','permute','true');
fdr = mafdr(pvalues,'BHFDR','true');
wtL_comp(:,4) = pvalues; 
wtL_comp(:,5) = fdr;

allbootL = wtL_comp(:,3);
allbootLP = wtL_comp(:,4);
meanDiffL = wtL_comp(:,1)-wtL_comp(:,2);

imagesc(meanDiffL'); colormap(cmapdiff); caxis([-.02 .02]);
for b = 1:162
    if allbootLP(b)>0 & allbootLP(b)<.05 & abs(meanDiffL(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,1,20,'k','filled')
    end
end
for b = 1:162
    if allbootL(b)>0 & allbootL(b)<.05 & abs(meanDiffL(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,1.25,20,'k','filled')
    end
end
for b = 1:162
    if allbootL(b)>0 & allbootL(b)<.05 & allbootLP(b)>0 & allbootLP(b)<.05 & abs(meanDiffL(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,.75,20,'k','filled')
    end
end









nsamps = 10000;
mwtS = mean(allwtS,'omitnan');
mamphS = mean(allamphS,'omitnan');
mpartS = mean(allpartS,'omitnan');

wtS_comp = zeros(162,5);
for b = 1:162
    [m1,m2,xt] = findBootMean(allwtS(:,b),allamphS(:,b),nsamps);
    wtS_comp(b,1) = m1; wtS_comp(b,2) = m2; wtS_comp(b,3) = xt./nsamps;
end
pvalues = mattest(allwtS',allamphS','permute','true');
fdr = mafdr(pvalues,'BHFDR','true');
wtS_comp(:,4) = pvalues; 
wtS_comp(:,5) = fdr;

allbootS = wtS_comp(:,3);
allbootSP = wtS_comp(:,4);
meanDiffS = wtS_comp(:,2)-wtS_comp(:,1);
figure(2); 
imagesc(meanDiffS'); colormap(cmapdiff); caxis([-.02 .02]);
for b = 1:162
    if allbootSP(b)>0 & allbootSP(b)<.05 & abs(meanDiffS(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,1,20,'k','filled')
    end
end
for b = 1:162
    if allbootS(b)>0 & allbootS(b)<.05 & abs(meanDiffS(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,1.25,20,'k','filled')
    end
end
for b = 1:162
    if allbootS(b)>0 & allbootS(b)<.05 & allbootSP(b)>0 & allbootSP(b)<.05 & abs(meanDiffS(b))>(50/90000)
        %hold on; text(b-.5,group+.1,'*');
        hold on; scatter(b,.75,20,'k','filled')
    end
end










%% scatter/swarm
y1 = allwtL(:,100);
y2 = allamphL(:,100);

x1 = ones(1,24);
x2 = 2+ones(1,6);

s1 = swarmchart(x1,y1,50,[.5 .5 .5],'filled'); hold on
s2 = swarmchart(x2,y2,50,[.5 .5 .5],'filled');
s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];

%%






nsamps = 100000;
GS = cell(5,1); GC = cell(5,1); GL = cell(5,1);
for group = 1:5
    group
tic
x0l = find(ratGROUP==group & isSOC==0 & ratGEN==0); [x0ls, x0lsi] = sort(ratGID(x0l));
x1l = find(ratGROUP==group & isSOC==0 & ratGEN==1); [x1ls, x1lsi] = sort(ratGID(x1l));
x0s = find(ratGROUP==group & isSOC==1 & ratGEN==0 & ratPGEN==0); [x0ss, x0ssi] = sort(ratGID(x0s));
xwt1s = find(ratGROUP==group & isSOC==1 & ratGEN==0 & ratPGEN==1); [xwt1ss, xwt1si] = sort(ratGID(xwt1s));
xko1s = find(ratGROUP==group & isSOC==1 & ratGEN==1 & ratPGEN==0); [xko1ss, xko1si] = sort(ratGID(xko1s));
x2s = find(ratGROUP==group & isSOC==1 & ratGEN==1 & ratPGEN==1); [x2ss, x2ssi] = sort(ratGID(x2s));


allwtL = ratCLASSH_F(x0l,orderedCZCoarse);
allkoL = ratCLASSH_F(x1l,orderedCZCoarse);
if group==1
   allwtL(13,:) = [];
end
wtkoL_comp = zeros(163,5);
for b = 1:163
    [m1,m2,xt] = findBootMean(allwtL(:,b),allkoL(:,b),nsamps);
    wtkoL_comp(b,1) = m1; wtkoL_comp(b,2) = m2; wtkoL_comp(b,3) = xt./nsamps;
end
pvalues = mattest(allwtL',allkoL','permute','true');
fdr = mafdr(pvalues,'BHFDR','true');
wtkoL_comp(:,4) = pvalues; wtkoL_comp(:,5) = fdr;
GL{group,1} = wtkoL_comp;



allwtSOC = ratSOCHIST_Z(x0s,socOrderZ);
allkoSOC = ratSOCHIST_Z(x2s,socOrderZ);
wtkoSOC_comp = zeros(161,5);
for b = 1:161
    [m1, m2, xt] = findBootMean(allwtSOC(:,b),allkoSOC(:,b),nsamps);
    wtkoSOC_comp(b,1) = m1; wtkoSOC_comp(b,2) = m2; wtkoSOC_comp(b,3) = xt./nsamps;
end
pvalues = mattest(allwtSOC',allkoSOC','permute','true');
fdr = mafdr(pvalues,'BHFDR','true');
wtkoSOC_comp(:,4) = pvalues; wtkoSOC_comp(:,5) = fdr;
GS{group,1} = wtkoSOC_comp;

allwtC = ratCLASSH_F(x0s,orderedCZCoarse);
allkoC = ratCLASSH_F(x2s,orderedCZCoarse);
wtkoC_comp = zeros(163,5);
for b = 1:163
    [m1, m2, xt] = findBootMean(allwtC(:,b),allkoC(:,b),nsamps);
    wtkoC_comp(b,1) = m1; wtkoC_comp(b,2) = m2; wtkoC_comp(b,3) = xt./nsamps;
end
pvalues = mattest(allwtC',allkoC','permute','true');
fdr = mafdr(pvalues,'BHFDR','true');
wtkoC_comp(:,4) = pvalues; wtkoC_comp(:,5) = fdr;
GC{group,1} = wtkoC_comp;
toc
end





%% pca plots
% BLUE DENSITY PLOT
figure(2); imagesc(D); colormap(cmapB./255); caxis([0 3.5E-4]); axis equal off;
hold on; scatter(llbwb(:,2),llbwb(:,1),.5,'k')
set(gcf,'Position',[166 81 1265 938]);


% LONE/SOCIAL RATIO
LLBound = cell(163,1);
for b = 1:163
    bi = LL2==b;
    fxmask = imclose(bi,ones(3,3));
    bd = bwboundaries(fxmask);
    try
    LLBound{b} = bd{1};
    catch
    end
end

imagesc(LLratio); axis equal off; colormap(gray); hold on;
set(gcf,'Position',[166 81 1265 938]);
for b = 1:163
    try
    k = behLabels(b);
    scatter(LLBound{b}(:,2),LLBound{b}(:,1),1,colorsCoarse(k+1,:));
    catch
    end
end

imagesc(LLratio); axis equal off; colormap(gray); hold on;
set(gcf,'Position',[166 81 1265 938]);
for b = 1:163
    try
    k = labelsCZCoarse(b);
    scatter(LLBound{b}(:,2),LLBound{b}(:,1),1,colorsCoarse(k+1,:));
    catch
    end
end




%% Mutual Information
EMB_WT = cell(length(x0s),2);
CW_WT = zeros(length(x0s),90000);

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
    if rec+396<1003
        rec2 = rec+396;
    else
        rec2 = rec-396;
    end
    EMB_AMPH{i,1} = wCOARSE{rec};
    EMB_AMPH{i,2} = wCOARSE{rec2};
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
FC_AMPH(isnan(FC_AMPH)) = 0;

EMB_AMPH2 = cell(length(x2s),2);
CW_AMPH2 = zeros(length(x2s),90000);
for i = 1:length(x2s)
     rec = x2s(i);
    if rec+396<1003
        rec2 = rec+396;
    else
        rec2 = rec-396;
    end
    EMB_AMPH2{i,1} = wCOARSE{rec};
    EMB_AMPH2{i,2} = wCOARSE{rec2};
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
FC_AMPH2(isnan(FC_AMPH2)) = 0;

figure(1); imagesc(FC_WT); colormap(cmapdiff); caxis([-.07 .07]); axis equal off
saveas(gcf,'figs20230917/MI_WT.tif'); 
colorbar; axis on
saveas(gcf,'figs20230917/MI_WT_CB.tif'); 

figure(2); imagesc(FC_AMPH); colormap(cmapdiff); caxis([-.07 .07]); axis equal off
saveas(gcf,'figs20230917/MI_AMPH.tif'); 
colorbar; axis on
saveas(gcf,'figs20230917/MI_AMPH_CB.tif'); 

figure(3); imagesc(FC_AMPH2); colormap(cmapdiff); caxis([-.07 .07]); axis equal off
saveas(gcf,'figs20230917/MI_AMPH2.tif'); 
colorbar; axis on
saveas(gcf,'figs20230917/MI_AMPH2_CB.tif'); 



%% individual embedding for LONGEVANS (wt/amph/wt-soc/amph-soc/part-soc)


allwtL = histIndivF(x0l,behOrderZ);
allamphL = histIndivF(x1l,behOrderZ);
allpartL = histIndivF(x2l,behOrderZ);

allwtS = histIndivF(x0s,behOrderZ);
allamphS = histIndivF(x1s,behOrderZ);
allpartS = histIndivF(x2s,behOrderZ);

symmwtS = histSOCC(x0s,:);
symmamphS = histSOCC(x1s,:);
symmpartS = histSOCC(x2s,:);


for i =1:4
    figure(i)
ss0s = symmwtS(:,i);
ss1s = symmamphS(:,i);
ss2s = symmpartS(:,i);

x1 = ones(1,length(x0s));
x2 = 2+ones(1,length(x1s));
x3 = 4+ones(1,length(x2s));

s1 = swarmchart(x1,ss0s,50,[.3 .3 .3],'filled'); hold on
s2 = swarmchart(x2,ss1s,50,[.3 .3 .3],'filled');
%s3 = swarmchart(x3,ss2s,50,[.3 .3 .3],'filled');

s1.MarkerEdgeColor = [0 0 0];
s2.MarkerEdgeColor = [0 0 0];
axis([0 4 0 1])
%s3.MarkerEdgeColor = [0 0 0];
set(gcf,'Position',[432 578 271 398]);
axis off
%saveas(gcf,['figs20230910/wt_amph_scatters_' num2str(i) '.tiff'])
end



s4.MarkerEdgeColor = [0 0 0];
s5.MarkerEdgeColor = [0 0 0];





idse = [x0l; x1l; x0s; x1s; x2s];
allseF = histIndivF(idse,behOrderZ);
allseC = histIndivC(idse,:);

 [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(allseF);
 [COEFFc, SCOREc, LATENTc, TSQUAREDc, EXPLAINEDc] = pca(allseC);

xcolors = [repmat([0 0 0],[length(x0l) 1]);repmat([0 1 1],[length(x1l) 1]);...
    repmat([0 0 1],[length(x0s) 1]);repmat([1 0 0],[length(x1s) 1]);repmat([0 1 0],[length(x2s) 1])];
 
 h = scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),75,xcolors,'filled');

h = scatter3(SCOREc(:,1),SCOREc(:,2),SCOREc(:,3),75,xcolors,'filled');
h = scatter3(SCOREc(:,1),SCOREc(:,2),SCOREc(:,3),75,xcolors);

alpha = 1;
set(h, 'MarkerFaceAlpha', alpha)
set(h, 'MarkerEdgeColor','k')

figure(1);
b = bar(COEFFc(:,1));
b.FaceColor = 'flat';
for i = 1:9
    b.CData(i,:) = colorsCoarse(i+1,:);
end
axis([0 10 -1 1]); axis equal

figure(2)
b2 = bar(COEFFc(:,2));
b2.FaceColor = 'flat';
for i = 1:9
    b2.CData(i,:) = colorsCoarse(i+1,:);
end
axis([0 10 -1 1]); axis equal

figure(3)
b3 = bar(COEFFc(:,3));
b3.FaceColor = 'flat';
for i = 1:9
    b3.CData(i,:) = colorsCoarse(i+1,:);
end
axis([0 10 -1 1]); axis equal

















