function [pvalue, testStatistic, lme] = testGroupLme(C1, C2, aC1, aC2)
C1n = C1+1e-6; C2n = C2+1e-6;
c1 = mean(C1n); % this is WT
c2 = mean(C2n); % this is KO

%Call = [C1; C2];
mjs = JSDiv(c1,c2);

a1n = unique(aC1);
a2n = unique(aC2);
a1 = length(a1n);
a2 = length(a2n);

wtid = aC1;
koid = aC2;

wtb = C1;
kob = C2;

cond = [zeros(size(wtb,1),1); ones(size(kob,1),1)];
aid = [aC1; aC2];

val = [];
for aID=1:size(wtb,1)
    val = [val; JSDiv(wtb(aID,:),c1)];
end

for aID=1:size(kob,1)
    val = [val; JSDiv(kob(aID,:),c1)];
end


tbl = table(cond,aid,val,'VariableNames',{'Condition','Animal','Expression'});

lme = fitlme(tbl,'Expression~Condition+(1|Animal)');

testStatistic = mjs;
pvalue = lme.Coefficients.pValue(2);
