%%

load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/fileOrder.mat')
load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/results202302112.mat')
load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/results202302112info.mat')
load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/comboClustering.mat')

FC = cell(size(FILECOMBO));
for i = 1:length(FILECOMBO)
    FC{i} = FILECOMBO{i}(110:end);
end

addpath(genpath('/Users/ugne/Dropbox/MultiFlyAnalysis_OLD/tSNE'))
skeleton = load('/Users/ugne/Dropbox/Label4D/skeletons/rat23.mat');

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
skeleton.mcolor = scM;
% skeleton.joint_names = joint_names;


Axl = find(ratGROUP==1 & isSOC==0);
Ax0s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==0);
Ax2s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==1);

Sxl = find(ratGROUP==5 & isSOC==0);
Sx0s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==0);
Sx2s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==1);


allLx = [Axl; Sxl];
allSx = [Ax0s; Ax2s; Sx0s; Sx2s];
allS2x = [Ax0s(19:36); Ax0s(1:18); Ax2s(20:38); Ax2s(1:19); Sx0s(10:18); Sx0s(1:9); Sx2s(10:18); Sx2s(1:9)];

allFLx = FC(allLx);
allFSx = FC(allSx);
allFSx2 = FC(allS2x);

allL = cell(size(allLx));
allS = cell(size([allSx allS2x]));

for i = 1:length(allL)
    try
        i
        p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' allFLx{i}],'pred');
        pred = p1.pred;
        for x = 1:23
            for y = 1:3
                pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
            end
        end
        allL{i} = pred;
    catch
    end
end

for i = 1:length(allS)
    i
    try
        p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' allFSx{i}],'pred');
        pred = p1.pred;
        for x = 1:23
            for y = 1:3
                pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
            end
        end
        p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' allFSx2{i}],'pred');
        pred2 = p2.pred;
        for x = 1:23
            for y = 1:3
                pred2(:,y,x) = smooth(medfilt1(pred2(:,y,x),3),3);
            end
        end
        allS{i,1} = pred;
        allS{i,2} = pred2;
    catch
    end
end

wrxL = [allLx; allSx; allS2x];
wrL = wFINE(wrxL);

wrxS = allSx;
wrS = wsocFINE(wrxS);

allP = [allL(:); allS(:)];
MA = cell(size(allP));
for i = 1:length(allP)
    try
    i
    MA{i} = alignDannceNF(allP{i});
    catch
    end
end


%% animal-autonomous behaviors

%G = makeGroupsAndSegments(wrL,163,ones(290,1),25);
G1 = makeGroupsAndSegments(ratCZ(wrxL),163,ones(290,1),25);
G2 = makeGroupsAndSegments(ratCZ(wrxL),163,ones(290,1),20);
G3 = makeGroupsAndSegments(ratCZ(wrxL),163,ones(290,1),15);
G4 = makeGroupsAndSegments(ratCZ(wrxL),163,ones(290,1),10);
G5 = makeGroupsAndSegments(ratCZ(wrxL),163,ones(290,1),5);
G = cell(size(G5));
for i = 1:163
    if size(G1{i},1) > 16
        G{i} = G1{i};
    elseif size(G2{i},1) > 16
        G{i} = G2{i};
    elseif size(G3{i},1) > 16
        G{i} = G3{i};
    elseif size(G4{i},1) > 16
        G{i} = G4{i};
    else G{i} = G5{i};
    end
end

offx = 250;
offy = 0;
offz = 250;
coffx = zeros(4,4); coffz = zeros(4,4);
for i = 1:4
    for j = 1:4
        coffx(i,j) = (i-1)*offx;
        coffz(i,j) = (j-1)*offz;
    end
end
coffz = coffz;

for tm = 1:163 %163
    try
        blabel = find(behOrderZ==tm);
        btag = behOrderTagZ(blabel);
        g = G{tm};
        lgs = size(g,1);
        if lgs > 16
            newG = g(randperm(lgs,16),:);
        else
            newG = G;
        end
        idxx = [1 2 5 6 9 11 12 15 18 21 25 28 29 34 35 38]
        
        newG = g(idxx,:);


        lenG = 60;
        gMarkers = cell(1,size(newG,1));

        for i = 1:size(newG,1)
            repG = MA{newG(i,1)}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers{i} = tgm(1:lenG,:,:);
            end
        end
        allM = []; allJ = []; allC = []; allmc = [];
        for i = 1:size(newG,1)
            NM = gMarkers{i};
            NM(:,1,:) = NM(:,1,:)+coffx(i);
            NM(:,3,:) = NM(:,3,:)+coffz(i);
            allM = cat(3,allM,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ = cat(1,allJ,newJ);
            allC = cat(1,allC,zeros(23,3));
            allmc = cat(1,allmc,sc);
        end
        sk2.color = allC; sk2.joints_idx = allJ; sk2.mcolor = allmc;

        close all;
        findic = figure('Name','Rat Test');
        h = cell(1,1);
        h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'MarkerSize',30,'LineWidth',1.5);
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition', [0.2651 0.0768 0.4185 0.7339]);
        view(h{1}.getAxes,20,30); axis equal off
        axis([-300 1000 -200 200 -200 1000]);
        set(gca,'Color','w')
        set(gcf,'Color','white');
        
        savePath = ['/Users/ugne/Dropbox/A_SDANNCE/b_LLAC/' num2str(btag) '_' num2str(blabel) '_' num2str(tm) '_s.avi'];
        frames = 2:lenG;

        % Uncomment to write the Animation to video.
        h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
    catch
    end
end


%% joint behaviors

GJ1 = makeGroupsAndSegments(ratSOCEZ(wrxS),161,ones(110,1),100);
GJ2 = makeGroupsAndSegments(ratSOCEZ(wrxS),161,ones(110,1),75);
GJ3 = makeGroupsAndSegments(ratSOCEZ(wrxS),161,ones(110,1),50);
GJ4 = makeGroupsAndSegments(ratSOCEZ(wrxS),161,ones(110,1),25);
GJ5 = makeGroupsAndSegments(ratSOCEZ(wrxS),161,ones(110,1),15);
GJ = cell(size(GJ5));

for i = 1:161
    if size(GJ1{i},1) > 16
        GJ{i} = GJ1{i};
    elseif size(GJ2{i},1) > 16
        GJ{i} = GJ2{i};
    elseif size(GJ3{i},1) > 16
        GJ{i} = GJ3{i};
    elseif size(GJ4{i},1) > 16
        GJ{i} = GJ4{i};
    else GJ{i} = GJ5{i};
    end
end


offx = 800;
offy = 800;
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


for tm = 1:161
    try
        blabel = find(socOrderZ==tm);
        btag = socOrderTagZ(blabel);
        G = GJ{tm};
        lgs = size(G,1);
        if lgs>9
            newG = G(randperm(lgs,9),:);
        else
            newG = G;
        end
        txx = [9 15 19 22 24 25 28 32 36];
        %txx = 28:36;    
        newG = G(txx,:);
        lenG = 150;
        gMarkers1 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = allS{newG(i,1),1}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers1{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers1{i} = tgm(1:lenG,:,:);
            end
        end
        gMarkers2 = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = allS{newG(i,1),2}(newG(i,2):newG(i,3),:,:);
            if size(repG,1) > lenG
                gMarkers2{i} = repG(1:lenG,:,:);
            else
                tgm = repmat(repG,[ceil(lenG/size(repG,1)) 1 1]);
                gMarkers2{i} = tgm(1:lenG,:,:);
            end
        end
        
        %transform gMarkers2 to centered coord
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
        h{1} = Keypoint3DAnimator(allM,sk2,'lineWidth',1.5);
        h{2} = Keypoint3DAnimator(allM2,sk3,'Axes',h{1}.Axes,'lineWidth',1.5);
        
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition',[0.3934 0.0679 0.5402 0.8232]);
                view(h{1}.getAxes,20,30); axis equal
                

        x1 = 0; y1 = 0;
        for i = 1:3
            x = x1+coffx(i,1);
            for j = 1:3
                y = y1+coffy(1,j);
                r  = 350;
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + x;
                yunit = r * sin(th) + y;
                h2 = plot(xunit, yunit,'Color',[0 0 0]); hold on;
            end
        end
         set(gcf,'Units','Normalized','OuterPosition',[0.3934 0.0679 0.5402 0.8232]);
                view(h{1}.getAxes,20,30); axis off equal
        % axis([-200 2200 -200 2200 -200 300]);
        axis([-400 2000 -400 2000 -200 300]);
        set(gcf,'Color','white');

        savePath = ['/Users/ugne/Dropbox/A_SDANNCE/b_LLJC4/' num2str(btag) '_' num2str(blabel) '_' num2str(tm) '_s.avi'];
        frames = 1:lenG;

        % Uncomment to write the Animation to video.
        h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
        pause(2)
    end
end


%% aligned -social perspective

MAS = cell(size(allS));
for i = 1:length(allS)
    try
    i
    [MAS{i,1} MAS{i,2}] = alignDannceSOC(allS{i,1},allS{i,2});
    catch
    end
end


[m1 m2] = alignDannceSOC(allS{i,1},allS{i,2});
h = cell(1,2);
h{1} = Keypoint3DAnimator(m1,skeleton,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(m2,skeleton,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
Animator.linkAll(h)
set(gca,'Color','k')

axis equal;
axis([-350 450 -350 450 -100 300]);
set(gcf,'Color','black');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);




