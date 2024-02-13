function scoreVal = scoreBetaHazardSimple(parasLower,xall,zall,p)
% groupIndex: 1: transplant, 0: nontransplant
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
scoreVal = length(parasLower(:));
% read bandwidth and bandwidth adjustment
% calculate \bb\trans\x
    btransx = xall*parasLower; % paras1 include upper d-by-d block.
    % calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
    [temp1,~,~,temp2] = trueMRL(zall,btransx,'mimicN');
    lambdaVec = (1./temp1-temp1).*temp2;

% lambdaVec(isnan(lambdaVec))=0;
lambdaVec(isinf(lambdaVec))=0;
    for j = 1:p
        scoreVal(j) = sum(lambdaVec.*xall(:,j));
    end
% scoreVal = scoreVal/samplesize;

end