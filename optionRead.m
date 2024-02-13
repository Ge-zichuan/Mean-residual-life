function [filefieldname,filefieldval] = optionRead(filenum)
% read file
filestr = fileread(strcat('options',num2str(filenum),'.txt'));
% break it into lines
filebyline = regexp(filestr, '\n', 'split');
filebyline = strtrim(filebyline);
% remove empty lines
filebyline(cellfun(@isempty,filebyline)) = [];
% split by ':'
% filefield = regexp(filebyline, ':', 'split');
filefieldname = {'p','d','samplesize','Dist of X','parameter','Dist of NonTransplant','Dist of waiting','Dist of Transplant','Dist of censoring','Censoring Rate'};
% filefieldname = [];
filefieldval = filebyline;
% for i = 1:length(filefieldname)
% %     filefieldname = [filefieldname{i};string(filefield{i}{1,1})];
%     filefieldval = [filefieldval;filebyline{i}];
% end
end