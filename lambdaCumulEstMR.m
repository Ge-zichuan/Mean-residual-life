function lambdaCumul = lambdaCumulEstMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,d,group)
estSample = length(inputg(:,1));
lambdaCumul = zeros(estSample,1);
uniqueT = unique(inputt);
for i = 1:length(uniqueT)
    lambdaCumul(inputt==uniqueT(i),1) = integral(@(x)lambdaEstMRSpec(parasLower,x,inputg(inputt==uniqueT(i),:),xall,zall,deltaall,samplesize,bands,d,group)...
        ,0,uniqueT(i),'ArrayValued',true);
end

end