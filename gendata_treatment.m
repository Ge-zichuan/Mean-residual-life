function [aout,probability] = gendata_treatment(factorF,coefficientA)
% factorF is the factor matrix of F, n-by-q
% coefficientA is the coefficient of F in the factor analysis (logistic
% regression), q-by-1.
probability = logsig(factorF*coefficientA);
% probability = probability*0.8+0.1;
aout = binornd(1,probability);
end