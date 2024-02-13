function lambdaDouInt = lambdaDouIntEstMR(parasLower,inputt1,inputg1,xall,zall,deltaall,samplesize,bands,d,group)
    innerCumulHazard = integral(@(x)lambdaEstMRSpec(parasLower,x,inputg1,xall,zall,deltaall,samplesize,bands,d,group)...
        ,0,inputt1,'ArrayValued',true);
    lambdaDouInt = exp(-innerCumulHazard);
end