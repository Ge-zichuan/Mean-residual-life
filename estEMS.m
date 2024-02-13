function expectMeanSurvival = estEMS(lambda,sortedZI,timeInterval,goodI)
% estimation of expected mean survival
sortedlambdaTilde = lambda(:,sortedZI);
lambdaCumTilde = cumsum(sortedlambdaTilde,2);
expectMeanSurvival = sum(timeInterval(:,1:length(goodI)).*exp(-lambdaCumTilde(:,1:length(goodI))),2);
end