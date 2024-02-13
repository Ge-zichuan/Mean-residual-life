function vic = calvicMR(paras0,xall,zdiffallN,zdiffallT,zdiffall,...
    kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p,d,bands)
% calculate vic
% Refer to: Yanyuan Ma, Xinyu Zhang,
% A validated information criterion to determine the structural dimension in dimension reduction models
vec1 = ones(p-1-d);
vec2 = zeros(p-1-d);
paraExp1 = zeros(p-d-1,d+1);
paraExp2 = zeros(p-d-1,d+1);
for i = 1:(p-d-1)
    for j = 1:d
        paraExp1(i,j) = paras0(i+1,j)-vec1(i)*paras0(1,j);
        paraExp2(i,j) = paras0(i+1,j)-vec2(i)*paras0(1,j);
    end
    paraExp1(i,d+1) = vec1(i);
    paraExp2(i,d+1) = vec2(i);
end

scoreValnew1 = scoreBetaMR(paraExp1,xall,zdiffallN,zdiffallT,zdiffall,...
    kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p-1,d+1,bands);
scoreValnew2 = scoreBetaMR(paraExp2,xall,zdiffallN,zdiffallT,zdiffall,...
    kbzN,kbzT,kbzND,kbzTD,deltaallN,deltaallT,samplesize,samplesizeN,samplesizeT,groupIndex,p-1,d+1,bands);
vic = (sum(scoreValnew1.^2)+sum(scoreValnew2.^2))...
    *sqrt(samplesize)/2.0+p*d*log(samplesize);
end

