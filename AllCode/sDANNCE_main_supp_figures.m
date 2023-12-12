%% sDANNCE all data supp


strain = 1;
figure(1);
allASDlfold = [ASD_LFOLD{1,1}; zeros(2,162); ASD_LFOLD{2,1}; zeros(2,162); ASD_LFOLD{3,1}; zeros(2,162); ...
ASD_LFOLD{4,1}; zeros(2,162); ASD_LFOLD{5,1};]; 
imagesc([allASDlfold]); caxis([-5 5]); colormap(cmapdiff); axis equal off

figure(2);
allASDsfold = [ASD_LFOLD{1,2}; zeros(2,162); ASD_LFOLD{2,2}; zeros(2,162); ASD_LFOLD{3,2}; zeros(2,162); ...
ASD_LFOLD{4,2}; zeros(2,162); ASD_LFOLD{5,2};]; 
imagesc([allASDsfold]); caxis([-5 5]); colormap(cmapdiff); axis equal off

figure(3);
allASDsocfold = [ASD_LFOLD{1,3}; zeros(2,156); ASD_LFOLD{2,3}; zeros(2,156); ASD_LFOLD{3,3}; zeros(2,156); ...
ASD_LFOLD{4,3}; zeros(2,156); ASD_LFOLD{5,3};]; 
imagesc([allASDsocfold]); caxis([-5 5]); colormap(cmapdiff); axis equal off





































