function m = mEstTMR(parasLower,inputt,inputg,xall,zall,wall,deltaall,samplesize,bands,d,upperLimit,group)
estSample = length(inputg(:,1));
uniqueT = unique(inputt);
intgral1 = lambdaCumulEstMR(parasLower,inputt,inputg,xall,zall,deltaall,samplesize,bands,d);
intgral2 = zeros(estSample,1);
for i = 1:length(uniqueT)
    intgral2(inputt==uniqueT(i)) = integral(@(x)lambdaDouIntEstMR(parasLower,x,inputg(inputt==uniqueT(i),:),xall,zall,deltaall,samplesize,bands,d,group)...
        ,uniqueT(i),upperLimit,'ArrayValued',true);
end
m = exp(intgral1).*intgral2;
end