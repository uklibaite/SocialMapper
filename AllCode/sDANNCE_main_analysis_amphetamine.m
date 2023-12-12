% sDANNCE_main_analysis_amphetamine

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





h = hist(behOrderTagZ(nupAP),1:9); h = hist(behOrderTagZ(ndownAP),1:9);
h = hist(socOrderTagZ(nupSOCA),0:8); h([1 4 2 5 3 8 9])


[mdiffLA, mfoldLA, fdrmLA, idcLA, nupLA, ndownLA] = testBehaviorSet5(allwtL, allamphL ,nnum);
[mdiffLS, mfoldLS, fdrmLS, idcLS, nupLS, ndownLS] = testBehaviorSet5(allwtL, allwtS ,nnum);

[mdiffSA, mfoldSA, fdrmSA, idcSA, nupSA, ndownSA] = testBehaviorSet5(allwtS, allamphS ,nnum);
[mdiffSP, mfoldSP, fdrmSP, idcSP, nupSP, ndownSP] = testBehaviorSet5(allwtS, allpartS ,nnum);
[mdiffAP, mfoldAP, fdrmAP, idcAP, nupAP, ndownAP] = testBehaviorSet5(allpartS, allamphS ,nnum);

[mdiffSOCA, mfoldSOCA, fdrmSOCA, idcSOCA, nupSOCA, ndownSOCA] = testBehaviorSet5(allwtSOC, allamphSOC ,nnum);
[mdiffSOCP, mfoldSOCP, fdrmSOCP, idcSOCP, nupSOCP, ndownSOCP] = testBehaviorSet5(allwtSOC, allpartSOC ,nnum);
[mdiffSOCAP, mfoldSOCAP, fdrmSOCAP, idcSOCAP, nupSOCAP, ndownSOCAP] = testBehaviorSet5(allpartSOC, allamphSOC ,nnum);

imagesc(mfoldLA'); colormap(cmapdiff); caxis([-3 3]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:162
    if ismember(b,idcLA)
        hold on; scatter(b,1,20,'k','filled');
    end
end

imagesc(mfoldLS'); colormap(cmapdiff); caxis([-3 3]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:162
    if ismember(b,idcLS)
        hold on; scatter(b,1,20,'k','filled');
    end
end

imagesc(mfoldSA'); colormap(cmapdiff); caxis([-3 3]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:162
    if ismember(b,idcSA)
        hold on; scatter(b,1,20,'k','filled');
    end
end

imagesc(mfoldSP'); colormap(cmapdiff); caxis([-3 3]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:162
    if ismember(b,idcSP)
        hold on; scatter(b,1,20,'k','filled');
    end
end

imagesc(mfoldSOCA'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:156
    if ismember(b,idcSOCA)
        hold on; scatter(b,1,20,'k','filled');
    end
end

imagesc(mfoldSOCP'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [62 370 1410 36]); axis off
for b = 1:156
    if ismember(b,idcSOCP)
        hold on; scatter(b,1,20,'k','filled');
    end
end


% 
dWT = zeros(24,65);
dAMPH = zeros(21,65);
for i = 1:24
    dWT(i,:) = hist(alldcc{x0s(i)},0:10:640)./90000;
end
for i = 1:21
    dAMPH(i,:) = hist(alldcc{x1s(i)},0:10:640)./90000;
end

% LONE (
% per rat change (lone wt->amph) 
% per rat change (soc. wt -> soc. amphetamine)
% per rat change (soc. wt -> soc. partner)
ratIDN = zeros(size(ratID));
for i = 1:length(ratIDN)
    p = ratID{i};
    ratIDN(i) = str2num(p(2:end));
end

idx0l = ratIDN(x0l)+1; idx0s = ratIDN(x0s)+1; 
idx1l = ratIDN(x1l)+1; idx1s = ratIDN(x1s)+1; 
idx2l = ratIDN(x2l)+1; idx2s = ratIDN(x2s)+1; 

alldcc = cell(1002,1);
alldcc(211:end) = [allDC; allDC];
perRatDist = cell(6,3); % 6 rats, 3 conditions
for rn = 1:6
    rnx = rn-1;
    xi0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0 & ratIDN==rnx);
    xi1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0 & ratIDN==rnx);
    xi0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1 & ratIDN==rnx);
    xi1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1 & ratIDN==rnx);
    xi2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1 & ratIDN==rnx);
    pf0 = PROXFRAC(xi0s);
    pf1 = PROXFRAC(xi1s);
    pf2 = PROXFRAC(xi2s);
    perRatDist{rn,1} = pf0;
    perRatDist{rn,2} = pf1;
    perRatDist{rn,3} = pf2;
end


allwtL_diff = zeros(6,162);

% S -> amphetamine social, S2 -> social partner (amph)

ratwtL = zeros(6,162); ratamphL = zeros(6,162); ratdiffL = zeros(6,162); ratfoldL = zeros(6,162);
ratwtS = zeros(6,162); ratamphS = zeros(6,162); ratdiffS = zeros(6,162); ratfoldS = zeros(6,162);
ratwtS2 = zeros(6,162); ratamphS2 = zeros(6,162); ratdiffS2 = zeros(6,162); ratfoldS2 = zeros(6,162);

ratwtSOC = zeros(6,156); ratamphSOC = zeros(6,156); ratdiffSOC = zeros(6,156); ratfoldSOC = zeros(6,156);
ratwtSOC2 = zeros(6,156); ratamphSOC2 = zeros(6,156); ratdiffSOC2 = zeros(6,156); ratfoldSOC2 = zeros(6,156);
ratfoldLS = zeros(6,162);
% group changes

% individual changes
for rn = 1:6
    rnx = rn-1;
    xi0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0 & ratIDN==rnx);
    xi1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0 & ratIDN==rnx);
    xi0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1 & ratIDN==rnx);
    xi1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1 & ratIDN==rnx);
    xi2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1 & ratIDN==rnx);
    

    l0l = histIndivF(xi0l,behOrderZ)+1e-6;
    l1l = histIndivF(xi1l,behOrderZ)+1e-6;
    s0s = histIndivF(xi0s,behOrderZ)+1e-6;
    s1s = histIndivF(xi1s,behOrderZ)+1e-6;
    s2s = histIndivF(xi2s,behOrderZ)+1e-6;
    soc0s = histJointF(xi0s,socOrderZ)+1e-6;
    soc1s = histJointF(xi1s,socOrderZ)+1e-6;
    soc2s = histJointF(xi2s,socOrderZ)+1e-6;
    
    ratwtL(rn,:) = mean(l0l);
    ratamphL(rn,:) = l1l;
    ratdiffL(rn,:) = l1l-mean(l0l); ratfoldL(rn,:) = log2(l1l./mean(l0l));

    

    ratwtS(rn,:) = mean(s0s);
    ratfoldLS(rn,:) = log2(mean(s0s)./mean(l0l));

    
    ratamphS(rn,:) = mean(s1s);
    ratdiffS(rn,:) = mean(s1s)-mean(s0s); ratfoldS(rn,:) = log2(mean(s1s)./mean(s0s));
    ratwtS2(rn,:) = mean(s0s);
    ratamphS2(rn,:) = mean(s2s);
    ratdiffS2(rn,:) = mean(s2s)-mean(s0s); ratfoldS2(rn,:) = log2(mean(s2s)./mean(s0s));
    
    ratwtSOC(rn,:) = mean(soc0s);
    ratamphSOC(rn,:) = mean(soc1s);
    ratdiffSOC(rn,:) = mean(soc1s)-mean(soc0s); ratfoldSOC(rn,:) = log2(mean(soc1s)./mean(soc0s));

    ratwtSOC2(rn,:) = mean(soc0s);
    ratamphSOC2(rn,:) = mean(soc2s);
    ratdiffSOC2(rn,:) = mean(soc2s)-mean(soc0s); ratfoldSOC2(rn,:) = log2(mean(soc2s)./mean(soc0s));
    
end


imagesc([ratwtL; ratamphL]); colormap(cmap1); 



figure(1); imagesc(ratfoldL); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(2); 
imagesc(mfoldLA'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:162
    if ismember(b,idcLA)
        hold on; scatter(b,1,20,'k','filled');
    end
end


figure(3); imagesc(ratfoldLS); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(4);
imagesc(mfoldLS'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:162
    if ismember(b,idcLS)
        hold on; scatter(b,1,20,'k','filled');
    end
end

figure(5); imagesc(ratfoldS); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(6);
imagesc(mfoldSA'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:162
    if ismember(b,idcSA)
        hold on; scatter(b,1,20,'k','filled');
    end
end

figure(11); imagesc(ratfoldS2); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(12);
imagesc(mfoldSP'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:162
    if ismember(b,idcSP)
        hold on; scatter(b,1,20,'k','filled');
    end
end

figure(7); imagesc(ratfoldSOC); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(8); imagesc(mfoldSOCA'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:162
    if ismember(b,idcSOCA)
        hold on; scatter(b,1,20,'k','filled');
    end
end

figure(9); imagesc(ratfoldSOC2); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position',[270 820 1386 74]); axis off
figure(10);
imagesc(mfoldSOCP'); colormap(cmapdiff); caxis([-5 5]); set(gcf,'Position', [270 820 1386 33]); axis off
for b = 1:156
    if ismember(b,idcSOCP)
        hold on; scatter(b,1,20,'k','filled');
    end
end







