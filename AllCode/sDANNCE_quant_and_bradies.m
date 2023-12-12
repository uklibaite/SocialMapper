%% video - dataset sdannce


% figs_1012 -> c3_eth
% 

% amphetamine = 6 rats, (30 lone + 45 social) = 120
x0l = find(ratGROUP==6 & isAMPH==0 & isSOC==0);
x1l = find(ratGROUP==6 & isAMPH==1 & isSOC==0);

x0s = find(ratGROUP==6 & isAMPH==0 & isSOC==1);
x1s = find(ratGROUP==6 & isAMPH==1 & isSOC==1);
x2s = find(ratGROUP==6 & isAMPH==2 & isSOC==1);



% arid = 8 rats, (39 lone + 37) = 113
x0l = find(ratGROUP==1 & isSOC==0);
x0s = find(ratGROUP==1 & isSOC==1 & ratGEN==0 & ratPGEN==0);
x2s = find(ratGROUP==1 & isSOC==1 & ratGEN==1 & ratPGEN==1);


% chd8 = 8 rats, (40 lone + 36) = 112
x0l = find(ratGROUP==2 & isSOC==0);
x0s = find(ratGROUP==2 & isSOC==1 & ratGEN==0 & ratPGEN==0);
x2s = find(ratGROUP==2 & isSOC==1 & ratGEN==1 & ratPGEN==1);

% grinb = 6 rats, (30 lone + 18) = 66
x0l = find(ratGROUP==3 & isSOC==0);
x0s = find(ratGROUP==3 & isSOC==1 & ratGEN==0 & ratPGEN==0);
x2s = find(ratGROUP==3 & isSOC==1 & ratGEN==1 & ratPGEN==1);

% nrxn1 = 8 rats, (40 lone + 36) = 112
x0l = find(ratGROUP==4 & isSOC==0);
x0s = find(ratGROUP==4 & isSOC==1 & ratGEN==0 & ratPGEN==0);
x2s = find(ratGROUP==4 & isSOC==1 & ratGEN==1 & ratPGEN==1);

% scn2a = 6 rats, (30 lone + 18) = 66
x0l = find(ratGROUP==5 & isSOC==0);
x0s = find(ratGROUP==5 & isSOC==1 & ratGEN==0 & ratPGEN==0);
x2s = find(ratGROUP==5 & isSOC==1 & ratGEN==1 & ratPGEN==1);


useLLAC = 
useLLJC = 
usePred = 

M2 = cell(1002,1);
for i = 1:1002
    pred = predALL{i,1};
    M2{i} = alignDannceNF(pred);
end

%% remake bradies (correct wr)
[groups,~,~] = makeGroupsAndSegments(ratCZ(:),max(max(LL)),ones(1,1002),25);




for tm = 1:max(max(LL))
    try
        G = groups{tm};
        lgs = size(G,1);
        if lgs>20
            newG = G(randperm(lgs,20),:);
        else
            newG = G;
        end
        lenG = 60;
        gMarkers = cell(1,size(newG,1));
        for i = 1:size(newG,1)
            repG = predALL{newG(i,1)}(newG(i,2):newG(i,3),:,:);
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
            NM(:,1,:) = NM(:,1,:);
            NM(:,3,:) = NM(:,3,:);
            
            allM = cat(3,allM,NM);
            newJ = skeleton.joints_idx+(i-1)*23;
            allJ = cat(1,allJ,newJ);
            if newG(i,1)<130
            allC = cat(1,allC,.7*ones(23,3));
            else
                allC = cat(1,allC,ones(23,3));
            end
            allmc = cat(1,allmc,sc);
        end
        sk2.color = allC; sk2.joints_idx = allJ; sk2.mcolor = allmc;
        
        close all;
        findic = figure('Name','Rat Test');
        h = cell(1,1);
        h{1} = Keypoint3DAnimator(allM,sk2,'Position',[0 0 1 1],'MarkerSize',30,'LineWidth',1.5);
        Animator.linkAll(h);
        set(gcf,'Units','Normalized','OuterPosition', [0.2651 0.0768 0.4185 0.7339]);
        view(h{1}.getAxes,20,30); axis equal
        axis([-300 450 -300 500 -50 200]);
        set(gca,'Color','k')

        set(gcf,'Color','black');
        
        savePath = ['/home/ugne/Dropbox/A_SDANNCE/bigRun/b_MNEW/' num2str(tm) '_s.avi'];
        frames = 2:lenG;
        
        % Uncomment to write the Animation to video.
        h{1}.writeVideo(frames, savePath, 'FPS', 25, 'Quality', 70);
    catch
    end
end

