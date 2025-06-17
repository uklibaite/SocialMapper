function [lme_pvals, lme_fdrm] = findSigBeh(C1,C2,id1,id2)   
C1n = C1+1e-6; C2n = C2+1e-6;

cond = [zeros(size(C1n,1),1); ones(size(C2n,1),1)];
aid = [id1; id2];


    lme_pvals = [];
    for clustID=1:size(C1n,2)
        val = [C1n(:,clustID); C2n(:,clustID)];
        tbl = table(cond,aid,val,'VariableNames',{'Condition','Animal','Expression'});

        try
            lme = fitlme(tbl,'Expression~Condition+(1|Animal)');
            lme_pvals = [lme_pvals lme.Coefficients.pValue(2)];
        catch
            lme_pvals = [lme_pvals NaN];
        end
    end

    lme_fdrm = mafdr(lme_pvals,'BHFDR','true');

    display(['Multi-level model, number p < 0.05 *BEFORE* FDR correction:', num2str(length(find(lme_pvals<0.05)))])
    display(['Multi-level model, number p < 0.05 *AFTER* FDR correction:', num2str(length(find(lme_fdrm<0.05)))])
