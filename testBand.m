function [numOfBands,numOfParas] = testBand(filenum,bands,eflag,para,paraTrue,ibandpos)
% Simple code to record all bands adjustment process
display([para';paraTrue']);
display([eflag ibandpos bands(ibandpos)]);
numOfBands = length(bands(:));
numOfParas = length(para(:));
allRecord = [numOfBands numOfParas bands(:)' eflag para(:)' paraTrue(:)'];
recordFormat = [repmat('%f ',1,numOfBands+2*numOfParas+3) '\n'];
fid = fopen(strcat('testBand',num2str(filenum),'.txt'), 'a+');
fprintf(fid, recordFormat, allRecord);
fclose(fid);
end