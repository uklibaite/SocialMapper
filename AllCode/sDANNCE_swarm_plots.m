%% amph - usage per 7 social cats


%% per context mean proximity - or just per rat per context

allcdwt = combineCells(alldcc(x0s));
allcdamph = combineCells(alldcc(x1s));
checkpart = combineCells(alldcc(x2s));

hwt = hist(allcdwt,0:10:640);
hamph = hist(allcdamph,0:10:640);
hpart = hist(checkpart,0:10:640);

hnwt = hwt./sum(hwt); hnamph = hamph./sum(hamph);
cdwt = cumsum(hnwt);
cdamph = cumsum(hnamph);


for i = 1:24
    plot(dWT(i,:),'LineWidth',.25,'Color',[0 0 0]); hold on
end
for i = 1:21
    plot(dAMPH(i,:),'LineWidth',.25,'Color',[.7 .7 .7]);
end
plot(hnwt,'LineWidth',4,'Color',[0 0 0]); hold on; plot(hnamph,'LineWidth',4,'Color',[.5 .5 .5]);




[h,p,ks2stat] = kstest2(allcdwt,allcdamph);

C = perRatDist';
C = C(:);
cents = [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19 21 22 23];
cc = repmat([0 0 0; 1 0 0; 0 0 1],[6 1]);
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);


C = perRatDist(:,1:2)';
C = C(:);
cents = [1 2 4 5 7 8 10 11 13 14 16 17];
cc = repmat([0 0 0; .5 .5 .5],[6 1]);
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);



c1 = combineCells(perRatDist(:,1));
c2 = combineCells(perRatDist(:,2));



%% swarms for amph - lookup for 121 135 159
bi= find(socOrderZ==121);
c1 = allwtSOC(:,bi);
c2 = allamphSOC(:,bi);
C = [{c1} {c2}];
cents = [1 2];
cc = [0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);
mfoldSOCA(bi)

bi= find(socOrderZ==24);
c1 = allwtSOC(:,bi);
c2 = allamphSOC(:,bi);
C = [{c1} {c2}];
cents = [1 2];
cc = [0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);
mfoldSOCA(bi)

bi= find(socOrderZ==135);
c1 = allwtSOC(:,bi);
c2 = allamphSOC(:,bi);
C = [{c1} {c2}];
cents = [1 2];
cc = [0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);
mfoldSOCA(bi)

bi= find(socOrderZ==159);
c1 = allwtSOC(:,bi);
c2 = allamphSOC(:,bi);
C = [{c1} {c2}];
cents = [1 2];
cc = [0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);
mfoldSOCA(bi)




%% swarms for ASD - do change instead (per rat?)
b = 13;
c1 = WT_USAGE{1,3}(:,b); c2 = WT_USAGE{4,3}(:,b); c3 = WT_USAGE{5,3}(:,b);
b1 = ASD_USAGE{1,3}(:,b); b2 = ASD_USAGE{4,3}(:,b); b3 = ASD_USAGE{5,3}(:,b);

C = [{c1} {b1} {c2} {b2} {c3} {b3}];
cents = [1 2 4 5 7 8];
cc = [0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);



b = 66;
c1 = WT_USAGE{1,3}(:,b); c2 = WT_USAGE{4,3}(:,b); c3 = WT_USAGE{5,3}(:,b);
b1 = ASD_USAGE{1,3}(:,b); b2 = ASD_USAGE{4,3}(:,b); b3 = ASD_USAGE{5,3}(:,b);

C = [{c1} {b1} {c2} {b2} {c3} {b3}];
cents = [1 2 4 5 7 8];
cc = [0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);





b = 152;
c1 = WT_USAGE{1,3}(:,b); c2 = WT_USAGE{4,3}(:,b); c3 = WT_USAGE{5,3}(:,b);
b1 = ASD_USAGE{1,3}(:,b); b2 = ASD_USAGE{4,3}(:,b); b3 = ASD_USAGE{5,3}(:,b);

C = [{c1} {b1} {c2} {b2} {c3} {b3}];
cents = [1 2 4 5 7 8];
cc = [0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);


figure(2); 
b = 147;
c1 = WT_USAGE{1,3}(:,b); c2 = WT_USAGE{4,3}(:,b); c3 = WT_USAGE{5,3}(:,b);
b1 = ASD_USAGE{1,3}(:,b); b2 = ASD_USAGE{4,3}(:,b); b3 = ASD_USAGE{5,3}(:,b);

C = [{c1} {b1} {c2} {b2} {c3} {b3}];
cents = [1 2 4 5 7 8];
cc = [0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5];
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);



% all strains, just diff
b = 13;
C = cell(5,1);
for i = 1:5
    C{i} = ASD_DIFF{i,3}(:,b);
end
cents = 1:2:10;
cc = zeros(5,3);
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);


% all strains, just diff
b = 13;
C = cell(5,1);
for i = 1:5
    C{i} = ASD_FOLD{i,3}(:,b);
end
cents = 1:2:10;
cc = zeros(5,3);
figure(1); line([0 10],[0 0],'Color',[.5 .5 .5],'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);


b = 66;
C = cell(5,1);
for i = 1:5
    C{i} = ASD_FOLD{i,3}(:,b);
end
cents = 1:2:10;
cc = zeros(5,3);
figure(2); line([0 10],[0 0],'Color',[.5 .5 .5],'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);

b = 152;
C = cell(5,1);
for i = 1:5
    C{i} = ASD_FOLD{i,3}(:,b);
end
cents = 1:2:10;
cc = zeros(5,3);
figure(3); line([0 10],[0 0],'Color',[.5 .5 .5],'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);


b = 147;
C = cell(5,1);
for i = 1:5
    C{i} = ASD_FOLD{i,3}(:,b);
end
cents = 1:2:10;
cc = zeros(5,3);
figure(4); line([0 10],[0 0],'Color',[.5 .5 .5],'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.5,cc);

%
% all strains, just diff
figure(1);
b = 13;
C = cell(5,2);
for i = 1:5
    C{i,2} = ASD_FOLD{i,3}(:,b);
    C{i,1} = WT_FOLD{i,3}(:,b);
end
C = C'; C = C(:);
cents = [1 2 4 5 7 8 10 11 13 14];
%cc = zeros(5,3);
cc = [.5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; ];
line([0 15],[0 0],'Color',[.5 .5 .5],'LineWidth',2,'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.4,cc);
axis off

figure(2);
b = 66;
C = cell(5,2);
for i = 1:5
    C{i,2} = ASD_FOLD{i,3}(:,b);
    C{i,1} = WT_FOLD{i,3}(:,b);
end
C = C'; C = C(:);
cents = [1 2 4 5 7 8 10 11 13 14];
%cc = zeros(5,3);
cc = [.5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; ];
line([0 15],[0 0],'Color',[.5 .5 .5],'LineWidth',2,'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.4,cc);
axis off

figure(3);
b = 152;
C = cell(5,2);
for i = 1:5
    C{i,2} = ASD_FOLD{i,3}(:,b);
    C{i,1} = WT_FOLD{i,3}(:,b);
end
C = C'; C = C(:);
cents = [1 2 4 5 7 8 10 11 13 14];
%cc = zeros(5,3);
cc = [.5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; ];
line([0 15],[0 0],'Color',[.5 .5 .5],'LineWidth',2,'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.4,cc);
axis off

aridWT = WT_USAGE{1,3}(:,b);
aridKO = ASD_USAGE{1,3}(:,b);
scn2aWT = WT_USAGE{5,3}(:,b);
scn2aKO = ASD_USAGE{5,3}(:,b);
(mean(aridKO)-mean(aridWT))./mean(aridWT)
(mean(scn2aKO)-mean(scn2aWT))./mean(scn2aWT)

figure(4);
b = 147;
C = cell(5,2);
for i = 1:5
    C{i,2} = ASD_FOLD{i,3}(:,b);
    C{i,1} = WT_FOLD{i,3}(:,b);
end
C = C'; C = C(:);
cents = [1 2 4 5 7 8 10 11 13 14];
%cc = zeros(5,3);
cc = [.5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; .5 .5 .5; 0 0 0; ];
line([0 15],[0 0],'Color',[.5 .5 .5],'LineWidth',2,'LineStyle','--'); hold on
[outinfo,tt] = mySwarmNA(C,cents,.4,cc);
axis off



