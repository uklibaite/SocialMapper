%% sdannce movies
% 1 - good tracking movie
% 2 - 


load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/fileOrder.mat')
load('/Users/ugne/Dropbox/A_SDANNCE/updatedSaves/results202302112.mat')


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

skeleton.color = scM;%zeros(23,3); 
skeleton.joints_idx = joints_idx;
skeleton.mcolor = scM;
% skeleton.joint_names = joint_names;


xi = find(ratGROUP==6 & isSOC==0);

fLE_LONE = FC(xi);
fLE_LONE = reshape(fLE_LONE,[6 5]);
predLE_LONE = cell(6,5);
for rat = 1:6
    for day = 1:5
        day
        p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_LONE{rat,day}],'pred');
        pred = p1.pred;
        for x = 1:23
            for y = 1:3
                pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
            end
        end
        predLE_LONE{rat,day} = pred;
    end
end


close all
fIndic= figure('Name','Rat Test');
h = cell(1);
h{1} = Keypoint3DAnimator(predLE_LONE{1}, skeleton, 'Position', [0 0  1 1],'MarkerSize',30);
set(gca,'Color','k')
Animator.linkAll(h)
axis equal;
axis([-350 450 -350 450 -20 300]);
set(gcf,'Color','black');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);


offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6); 
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:6); coffx = coffx(1:5,1:6);
allM = []; allJ = []; allC = []; allmc = []; count = 0;
for rat= 1:6
    for day = 1:5
    cx = coffy(rat,day); cy = coffx(rat,day); cz = coffz(rat,day);
    NM = predLE_LONE{rat,day}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    
    newJ = skeleton.joints_idx+count*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
    count = count+1;
    end
end
sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
close all;
findic = figure('Name','Rat Test');
h = cell(1);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1);
axis equal off

unitcolor = cell(6,5);
for i = 1:6
    for j = 1:5
        if j~=4
            unitcolor{i,j} = [0 0 0];
        else
            unitcolor{i,j} = [0 1 1];
        end
    end
end

x1 = 45; y1 = 90;
for i = 1:6

    for j = 1:5
        x = x1+coffy(i,j);
        y = y1+coffx(i,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h2 = plot(xunit, yunit, 'Color',unitcolor{i,j},'LineWidth',2); hold on;
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 3595 -262 4440 -20 260]);
Animator.linkAll(h)
set(gcf,'Color','white');


savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allLE_LONE.avi'];
frames = 2:25:45000;

% Uncomment to write the Animation to video.
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);



%%

xi = find(ratGROUP==6 & isSOC==1 & isAMPH==0);

fLE_SOC1 = FC(xi(1:24));
fLE_SOC2 = FC(xi(25:end));


predLE_SOC = cell(24,2);
for rat = 1:24
    rat
    p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC1{rat}],'pred');
    pred = p1.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end

    predLE_SOC{rat,1} = pred;

    p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC2{rat}],'pred');
    pred = p2.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end
    predLE_SOC{rat,2} = pred;
end





p1 = predLE_SOC{1,1}; p2 = predLE_SOC{1,2};
h = cell(1,2);
h{1} = Keypoint3DAnimator(p1,skeleton,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(p2,skeleton,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
Animator.linkAll(h)
set(gca,'Color','k')

axis equal;
axis([-350 450 -350 450 -20 300]);
set(gcf,'Color','black');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);





offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6);
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:3); coffx = coffx(1:5,1:3);
allM = []; allJ = []; allC = []; allmc = [];
for i = 1:24
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,1}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
end

allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
for i = 1:24
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,2}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM2 = cat(3,allM2,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ2 = cat(1,allJ2,newJ);
    allC2 = cat(1,allC2,skeleton.color);
end

sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
sk2.mcolor = zeros(size(sk2.mcolor));
sk3 = sk2; sk3.mcolor = ones(size(sk2.mcolor)).*.7;
sk3.color = ones(size(sk2.color)).*.7;
close all;
findic = figure('Name','Rat Test');
h = cell(1,2);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(allM2,sk3,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
axis equal off
Animator.linkAll(h)

x1 = 90; y1 = 90;
for i = 1:6
    x = x1+coffx(i,1);
    for j = 1:4
        y = y1+coffy(1,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h2 = plot(xunit, yunit, 'Color',[0 0 1],'LineWidth',2); hold on;
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 4440 -262 2840 -25 260]);
Animator.linkAll(h)
set(gcf,'Color','white');
savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allWT_SOC.avi'];
frames = 2:25:90000;
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);


clear predLE_SOC
%%
xi1 = find(ratGROUP==6 & isSOC==1 & isAMPH==1);
%xi2 = find(ratGROUP==6 & isSOC==1 & isAMPH==2);
xi2 = zeros(size(xi1));
for i = 1:length(xi1)
    if xi1(i)-396<=211
        xi2(i) = xi1(i)+396;
    else
        xi2(i) = xi1(i)-396;
    end
end


fLE_SOC1 = FC(xi1);
fLE_SOC2 = FC(xi2);


predLE_SOCA = cell(21,2);
for rat = 1:21
    rat
    p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC1{rat}],'pred');
    pred = p1.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end

    predLE_SOCA{rat,1} = pred;

    p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC2{rat}],'pred');
    pred = p2.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end
    predLE_SOCA{rat,2} = pred;
end





p1 = predLE_SOCA{1,1}; p2 = predLE_SOCA{1,2};
h = cell(1,2);
h{1} = Keypoint3DAnimator(p1,skeleton,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(p2,skeleton,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
Animator.linkAll(h)
set(gca,'Color','k')

axis equal;
axis([-350 450 -350 450 -20 300]);
set(gcf,'Color','black');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);





offx = 800;
offy = 800;
offz = 0;
coffx = zeros(7,7); coffy = zeros(7,7); coffz = zeros(7,7);
for i = 1:7
    for j = 1:7
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:3); coffx = coffx(1:5,1:3);
allM = []; allJ = []; allC = []; allmc = [];
for i = 1:21
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOCA{i,1}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
end

allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
for i = 1:21
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOCA{i,2}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM2 = cat(3,allM2,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ2 = cat(1,allJ2,newJ);
    allC2 = cat(1,allC2,skeleton.color);
end

sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
sk2.mcolor = zeros(size(sk2.mcolor));
sk3 = sk2; sk3.mcolor = ones(size(sk2.mcolor)).*.7;
sk3.color = ones(size(sk2.color)).*.7;
close all;
findic = figure('Name','Rat Test');
h = cell(1,2);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(allM2,sk3,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
axis equal off
Animator.linkAll(h)

x1 = 90; y1 = 90;
for i = 1:7
    x = x1+coffx(i,1);
    for j = 1:3
        y = y1+coffy(1,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        if i == 5 && j==5

        else
        h2 = plot(xunit, yunit, 'Color',[1 0 0],'LineWidth',2); hold on;
        end
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-260 5240 -262 2040 -30 260]);
Animator.linkAll(h)
set(gcf,'Color','white');
savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allWT_AMPH.avi'];
frames = 2:25:90000;
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);

%%
xi1 = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==1);

fLE_SOC1 = FC(xi1(1:19));
fLE_SOC2 = FC(xi1(20:end));


predLE_SOC = cell(19,2);
for rat = 1:19
    rat
    p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC1{rat}],'pred');
    pred = p1.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end

    predLE_SOC{rat,1} = pred;

    p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC2{rat}],'pred');
    pred = p2.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end
    predLE_SOC{rat,2} = pred;
end

offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6);
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:3); coffx = coffx(1:5,1:3);
allM = []; allJ = []; allC = []; allmc = [];
for i = 1:18
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,1}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
end

allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
for i = 1:18
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,2}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM2 = cat(3,allM2,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ2 = cat(1,allJ2,newJ);
    allC2 = cat(1,allC2,skeleton.color);
end

sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
sk2.mcolor = zeros(size(sk2.mcolor));
sk3 = sk2; sk3.mcolor = ones(size(sk2.mcolor)).*.7;
sk3.color = ones(size(sk2.color)).*.7;
close all;
findic = figure('Name','Rat Test');
h = cell(1,2);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(allM2,sk3,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
axis equal off
Animator.linkAll(h)

x1 = 45; y1 = 90;
for i = 1:6
    x = x1+coffx(i,1);
    for j = 1:3
        y = y1+coffy(1,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h2 = plot(xunit, yunit, 'Color',[0 0 0],'LineWidth',2); hold on;
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 4440 -262 2040 -25 260]);
Animator.linkAll(h)
set(gcf,'Color','white');
savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allARIDKO_SOC.avi'];
frames = 2:25:90000;
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);

%%
xi1 = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==0);

fLE_SOC1 = FC(xi1(1:18));
fLE_SOC2 = FC(xi1(19:end));


predLE_SOC = cell(18,2);
for rat = 1:18
    rat
    p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC1{rat}],'pred');
    pred = p1.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end

    predLE_SOC{rat,1} = pred;

    p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC2{rat}],'pred');
    pred = p2.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end
    predLE_SOC{rat,2} = pred;
end

offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6);
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:3); coffx = coffx(1:5,1:3);
allM = []; allJ = []; allC = []; allmc = [];
for i = 1:18
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,1}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
end

allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
for i = 1:18
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{i,2}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM2 = cat(3,allM2,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ2 = cat(1,allJ2,newJ);
    allC2 = cat(1,allC2,skeleton.color);
end

sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
sk2.mcolor = zeros(size(sk2.mcolor));
sk3 = sk2; sk3.mcolor = ones(size(sk2.mcolor)).*.7;
sk3.color = ones(size(sk2.color)).*.7;
close all;
findic = figure('Name','Rat Test');
h = cell(1,2);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(allM2,sk3,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
axis equal off
Animator.linkAll(h)

x1 = 45; y1 = 90;
for i = 1:6
    x = x1+coffx(i,1);
    for j = 1:3
        y = y1+coffy(1,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h2 = plot(xunit, yunit, 'Color',[.7 .7 .7],'LineWidth',2); hold on;
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 4440 -262 2040 -25 260]);
Animator.linkAll(h)
set(gcf,'Color','white');
savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allARIDWT_SOC.avi'];
frames = 2:25:90000;
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);


%%


xi1 = find(ratGROUP==5 & isSOC==0);

fLE_LONE = FC(xi1);
fLE_LONE = reshape(fLE_LONE,[6 5]);
predLE_LONE = cell(6,5);
for rat = 1:6
    for day = 1:5
        day
        p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_LONE{rat,day}],'pred');
        pred = p1.pred;
        for x = 1:23
            for y = 1:3
                pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
            end
        end
        predLE_LONE{rat,day} = pred;
    end
end


close all
fIndic= figure('Name','Rat Test');
h = cell(1);
h{1} = Keypoint3DAnimator(predLE_LONE{1}, skeleton, 'Position', [0 0  1 1],'MarkerSize',30);
set(gca,'Color','k')
Animator.linkAll(h)
axis equal;
axis([-350 450 -350 450 -20 300]);
set(gcf,'Color','black');
% set(fIndic,'pos',pos);
view(h{1}.getAxes,-30,40);


offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6); 
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:6); coffx = coffx(1:5,1:6);
allM = []; allJ = []; allC = []; allmc = []; count = 0;
for rat= 1:6
    for day = 1:5
    cx = coffy(rat,day); cy = coffx(rat,day); cz = coffz(rat,day);
    NM = predLE_LONE{rat,day}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    
    newJ = skeleton.joints_idx+count*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
    count = count+1;
    end
end
sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
close all;
findic = figure('Name','Rat Test');
h = cell(1);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1);
axis equal off

unitcolor = cell(6,5);
for i = 1:6
    for j = 1:5
        if i <=3
            unitcolor{i,j} = [.7 .7 .7];
        else
            unitcolor{i,j} = [0 0 0];
        end
    end
end

x1 = 45; y1 = 90;
for i = 1:6

    for j = 1:5
        x = x1+coffy(i,j);
        y = y1+coffx(i,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h2 = plot(xunit, yunit, 'Color',unitcolor{i,j},'LineWidth',2); hold on;
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 3595 -262 4440 -20 260]);
Animator.linkAll(h)
set(gcf,'Color','white');


savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allSCN2A_LONE.avi'];
frames = 2:25:45000;

% Uncomment to write the Animation to video.
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);

%%
xi1 = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==1);
xi2 = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==0);

fLE_SOC1 = FC([xi1(1:9); xi2(1:9)]);
fLE_SOC2 = FC([xi1(10:end); xi2(10:end)]);


predLE_SOC = cell(18,2);
for rat = 1:18
    rat
    p1 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC1{rat}],'pred');
    pred = p1.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end

    predLE_SOC{rat,1} = pred;

    p2 = load(['/Volumes/olveczky_lab/Lab/dannce_rig/' fLE_SOC2{rat}],'pred');
    pred = p2.pred;
    for x = 1:23
        for y = 1:3
            pred(:,y,x) = smooth(medfilt1(pred(:,y,x),3),3);
        end
    end
    predLE_SOC{rat,2} = pred;
end

offx = 800;
offy = 800;
offz = 0;
coffx = zeros(6,6); coffy = zeros(6,6); coffz = zeros(6,6);
for i = 1:6
    for j = 1:6
        coffx(i,j) = (i-1)*offx;
        coffy(i,j) = (i-1)*offy;
        coffz(i,j) = (j-1)*offz;
    end
end   
coffy = coffy'; %coffy = coffy(1:5,1:3); coffx = coffx(1:5,1:3);
allM = []; allJ = []; allC = []; allmc = [];
order = [1 2 3 10 11 12 4 5 6 13 14 15 7 8 9 16 17 18];
for i = 1:18

    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{order(i),1}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM = cat(3,allM,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ = cat(1,allJ,newJ);
    allC = cat(1,allC,skeleton.color);
end

allM2 = []; allJ2 = []; allC2 = []; allmc2 = [];
for i = 1:18
    cx = coffx(i); cy = coffy(i); cz = coffz(i);
    NM = predLE_SOC{order(i),2}(1:90000,:,:);
    NM(:,1,:) = NM(:,1,:)+cx;
    NM(:,2,:) = NM(:,2,:)+cy;
    NM(:,3,:) = NM(:,3,:)+cz;
    allM2 = cat(3,allM2,NM);
    newJ = skeleton.joints_idx+(i-1)*23;
    allJ2 = cat(1,allJ2,newJ);
    allC2 = cat(1,allC2,skeleton.color);
end

sk2.mcolor = allC; sk2.joints_idx = allJ; sk2.color = zeros(size(allC));
sk2.mcolor = zeros(size(sk2.mcolor));
sk3 = sk2; sk3.mcolor = ones(size(sk2.mcolor)).*.7;
sk3.color = ones(size(sk2.color)).*.7;
close all;
findic = figure('Name','Rat Test');
h = cell(1,2);
h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'lineWidth',1.5);
h{2} = Keypoint3DAnimator(allM2,sk3,'Position',[0 0 1 1],'lineWidth',1.5,'Axes',h{1}.Axes);
axis equal off
Animator.linkAll(h)

x1 = 45; y1 = 90;
for i = 1:6
    x = x1+coffx(i,1);
    for j = 1:3
        y = y1+coffy(1,j);
        r  = 350;
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        if i <=3
        h2 = plot(xunit, yunit, 'Color',[0 0 0],'LineWidth',2); hold on;
        else
                   h2 = plot(xunit, yunit, 'Color',[.7 .7 .7],'LineWidth',2); hold on;
        end
    end
end
set(gcf,'Position',[296 62 1421 964]);
view(h{1}.getAxes,5,30); axis equal
axis([-305 4440 -262 2040 -25 260]);
Animator.linkAll(h)
set(gcf,'Color','white');
savePath = ['/Users/ugne/Dropbox/SDANNCE_MARCH_FIGS/allSCN2A_SOC.avi'];
frames = 2:25:90000;
h{1}.writeVideo(frames, savePath, 'FPS', 50, 'Quality', 70);


%% 
addpath(genpath('/Users/ugne/Dropbox/MultiFlyAnalysis_OLD/tSNE'))
load('/Users/ugne/Dropbox/Label4D/skeletons/rat23.mat')

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
sc1 = scM;
skeleton.color = scM;%zeros(23,3); 
skeleton.joints_idx = joints_idx;
skeleton.mcolor = scM;
skeleton.joint_names = joint_names;

projectFolderA = '/Users/ugne/Dropbox/A_SDANNCE/2022_09_22_M1_M2/'

pred1 = load([projectFolderA 'SDANNCE/bsl0.5_FM_rat1/save_data_AVG0.mat']);
pred2 = load([projectFolderA 'SDANNCE/bsl0.5_FM_rat2/save_data_AVG0.mat']);
p1 = pred1.pred;
p2 = pred2.pred;
pSOC = cat(3,p1,p2);

for i = 1:3
    for j = 1:46
        pSOC(:,i,j) = medfilt1(pSOC(:,i,j),3);
    end
end


calibPaths = collectCalibrationPaths(projectFolderA);
params = cellfun(@(X) {load(X)}, calibPaths);
vidName = '0.mp4';
vidPaths = collectVideoPaths(projectFolderA,vidName);
sync = cell(numel(calibPaths),1);
nFrames = 90000;
nKeypoints = 23;

for i = 1:numel(calibPaths)
    sync_struct = struct('data_2d',[],'data_3d',[],'data_frame',[],'data_sampleID',[]);
    sync_struct.data_2d = zeros(nFrames,nKeypoints*2);
    sync_struct.data_3d = zeros(nFrames,nKeypoints*3);
    sync_struct.data_frame = 0:(nFrames-1);
    sync_struct.data_sampleID = sync_struct.data_frame;
    sync{i} = sync_struct;
end

[c,orientations,locations] = deal(cell(6,1));
for i = 1:numel(c)
    % Get all parameters into cameraParameters object.
    cp = calibPaths{i}; load(cp);
    rotationVector = rotationMatrixToVector(r);
    translationVector =t;
    c{i} = cameraParameters('IntrinsicMatrix',K,...
        'ImageSize',[1200 1920],'RadialDistortion',RDistort,...
        'TangentialDistortion',TDistort,...
        'RotationVectors',rotationVector,...
        'TranslationVectors',translationVector);
    
    % Also save world location and orientation
    orientations{i} = r';
    locations{i} = -translationVector*orientations{i};
end


nCam = 1;
cpt = zeros(90000,23*2,2);
for f = 1:90000
    points = squeeze(pSOC(f,:,:))';
    iP = worldToImage(c{nCam}.Intrinsics,c{nCam}.RotationMatrices,c{nCam}.TranslationVectors,points,'ApplyDistortion',true);
    cpt(f,:,:) = iP;
end
cpt2 = permute(cpt,[1 3 2]);



ax = squeeze(pSOC(:,1,:)); ax = ax(:);
ay = squeeze(pSOC(:,2,:)); ay = ay(:);
az = squeeze(pSOC(:,3,:)); az = az(:);
minx = prctile(ax,.01); miny = prctile(ay,.01); minz = prctile(az,.01);
maxx = prctile(ax,99.9); maxy = prctile(ay,99.9); maxz = prctile(az,99.9);
limits = [minx-30 maxx+30 miny-30 maxy+30 minz-10 maxz+40];
minz = -5;
cyl.color = [.1 .5 .1];
cyl.h = 10;
cyl.r = 335;
cyl.n = 100;
cyl.cent = [40 95];
cyl.minz = minz;

frame = 1;
nVid = 1;
vr = VideoReader(vidPaths{nVid});
V = read(vr,[5001 8000]);

skelSOC.joints_idx = [skeleton.joints_idx; skeleton.joints_idx+23];
skelSOC.mcolor = [sc1; sc1];
skelSOC.color = zeros(46,3); skelSOC.color(1:23,:) = .5*ones(23,3);
skelSOC2 = skelSOC;
skelSOC2.color(1:23,:) = ones(23,3);
skelSOC2.color(24:end,:) = .7*ones(23,3);

skelSOC.color(1:23,:) = ones(23,3);
skelSOC.color(24:end,:) = .7*ones(23,3);
nVid = 1;
close all
h = cell(2,1);
% points = cpt2(starts(i):ends(i),:,:);
points = cpt2(5001:8000,:,:);
fMulti = figure('Name','Combination Test 2');
h{2} = Keypoint3DAnimatorObj(pSOC(5001:8000,:,:),skelSOC,cyl,limits,'Position',[.5 .25 .5 .5],'LineWidth',2,'MarkerSize',50);
set(gcf,'color','k');
axis equal off
view(locations{nVid})
h{1} = VideoAnimator(V*2.5,'Position',[0 .25 .5 .5]);
axis equal off
% Link all animators in the cell array.
Animator.linkAll(h)
set(gcf,'color','k');
axis equal off
set(gcf, 'Position',[586 89 1195 928]);

















