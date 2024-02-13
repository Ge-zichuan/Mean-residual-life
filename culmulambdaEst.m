function culmulambdaVec = culmulambdaEst(inputt,inputg,zall,aall,deltaall,samplesize,goodIndex,bands)
% This is the estimation of Lambda, integral of lambda from 0 to t.
culmulambdaVec = integral(@(x)lambdaEst(x,inputg,zall,aall,deltaall,samplesize,goodIndex,bands),0,inputt,'ArrayValued',true);
end