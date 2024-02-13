function [bandPlot] = testBandAnalysis(filenum)
% Simple code to record all bands adjustment process
allInformation = dlmread(strcat('testBand',num2str(filenum),'.txt'));
numberOfTrials = length(allInformation(:,1));
numOfBands = allInformation(1,1);
numOfParas = allInformation(1,2);
allBands = allInformation(:,3:(2+numOfBands));
allEflags = allInformation(:,2+numOfBands+1);
allParas = allInformation(:,(4+numOfBands):(3+numOfBands+numOfParas));
paraTrue = allInformation(1,(4+numOfBands+numOfParas):(3+numOfBands+2*numOfParas));
s = ' ';%n = '\n';
% for i = 1:20
%     if (length(unique(allBands(:,i))) ~= 1)
%         display(['Band',s,num2str(i),s,'has',s,num2str(length(unique(allBands(:,i)))),s,'values']);
%     end
% end
% whichBand = input('Which band?');
for whichBand = 1:27
% while (ismember(whichBand,1:20))
    otherBands = allBands(:,setdiff(1:numOfBands,whichBand));
    [C,ia,ic] = unique(otherBands,'rows');
    numOfGroups = length(unique(ic));
    for i = 1:numOfGroups
        rowIndex = (ic==i);
        bandChangeX = allBands(rowIndex,whichBand);
        paraChangeY = allParas(rowIndex,:);
        eflagChange = allEflags(rowIndex,:);
        %     subplot(1,numOfGroups,i);
        if (length(unique(bandChangeX))==1)
        else
            figure;
            hold on;
            [xo,xoi] = sort(bandChangeX);
            cla;
            colors = colormap(lines(8));
            position = 0;
            for j = 1:numOfParas
                %                 plot(xo,paraChangeY(xoi,j),'-','linewidth',2,'Color',colors(j,:)); hold on;
                %                 plot(xo,(ones(length(xo),1)*paraTrue(j)),'-.','linewidth',2,'Color',colors(j,:)); hold on;
                plot(xo,paraChangeY(xoi,j)-paraTrue(j),'-','linewidth',2,'Color',colors(j,:)); hold on;
                position = max(max(paraChangeY(xoi,j)-paraTrue(j)),position);
            end
            plot(xo,zeros(length(xo),1),'.','linewidth',2,'Color','k'); hold on;
            title(strcat('The',s,num2str(whichBand),' th',' band'));
            for j = 1:numel(xo)
                text(xo(j),position,num2str(eflagChange(j)));
            end
            titleBands = reshape(allBands(ia(i),:),4,[]);
            titleBands(whichBand) = 999;
            xlabel(num2str(titleBands));
            hold off;
%             s = waitforbuttonpress;
        end
    end
%     for i = 1:20
%         if (length(unique(allBands(:,i))) ~= 1)
%             display(['Band',s,num2str(i),s,'has',s,num2str(length(unique(allBands(:,i)))),s,'values']);
%         end
%     end
%     whichBand = input('Which band?');
% end
end
end