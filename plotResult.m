function plotT = plotResult(dataX,dataXZ,dataY,dataZ,baseLine,baseLineZ,optionsPara,displayType,groupType,saveName,bands)
switch groupType
    case 'N'
        plotT = figure('visible',displayType);%,'Position', get(0, 'Screensize')
        bbxallposs = unique(dataX(:,2:end),'rows');
        totalPlotN = length(bbxallposs);
        for i = 1:totalPlotN
            subplot(ceil(sqrt(totalPlotN)),ceil(sqrt(totalPlotN)),i)
            bbx = bbxallposs(i,:);
            plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows'));
            plot(dataX(plotsubindexT,1),baseLine(plotsubindexT),'-');hold on;
            scatter(dataX(plotsubindexT,1),mean(dataY(:,plotsubindexT),1),'filled');
            scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.975,1),'+');
            scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.025,1),'o');
            ylim([-1 max(10,1.2*max(baseLine(plotsubindexT)))]);
            if (i == 1)
                title(num2str(bands));hold off;
            else
                title(num2str(bbx));hold off;
            end
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_N',num2str(optionsPara(1)),'_',regexprep(num2str([optionsPara(2),optionsPara(3)]),'\s+','_'),'.png')],'png')
        %         close(plotT);
    case'T'
        plotT = figure('visible',displayType);%,'Position', get(0, 'Screensize')
        bbxallposs = unique(dataX(:,2:end),'rows');
        plotRows = 1:round(length(unique(bbxallposs(:,2)))/5):length(unique(bbxallposs(:,2)));
        plotCols = 1:round(length(unique(bbxallposs(:,1)))/5):length(unique(bbxallposs(:,1)));
        for i = 1:length(plotRows)
            for j = 1:length(plotCols)
                bbx = bbxallposs((plotCols(j)-1)*length(unique(bbxallposs(:,2)))+plotRows(i),:);
                plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows'));
                subplot(length(plotRows),length(plotCols),(i-1)*length(plotCols)+j);
                plot(dataX(plotsubindexT,1),baseLine(plotsubindexT),'-');hold on;
                scatter(dataX(plotsubindexT,1),nanmean(dataY(:,plotsubindexT),1),'filled');
                scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.975,1),'+');
                scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.025,1),'o');
                ylim([-1 max(10,5*max(baseLine(plotsubindexT)))]);
                if (i == 1 && j==1)
                    title(num2str(bands));hold off;
                else
                    title(num2str(bbx));hold off;
                end
                %     title(num2str([bands(1,:);bands(2,:)]));hold off;
            end
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_T',num2str(optionsPara(1)),'_',regexprep(num2str([optionsPara(2),optionsPara(3)]),'\s+','_'),'.png')],'png')
        %         close(plotT);
    case'improve'
        plotT = figure('visible',displayType,'Position', get(0, 'Screensize'));
        bbxallposs = unique(dataX(:,2:end),'rows');
        estW = unique(bbxallposs(:,end));
        estbbx = unique(bbxallposs(:,1:(end-1)));
        for i = 1:(length(estW)/2)
            for j = 1:(length(estbbx)/2)
                bbx = bbxallposs((j-1)*length(estW)+i,:);
                plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows'));
                yT = dataY(:,plotsubindexT);
                plotsubindexN = find(ismember(dataXZ(:,2:end), bbx(:,1:(end-1)),'rows'));
                yN = dataZ(:,plotsubindexN);
                ydiff = yT(:,1:(end-i))-yN(:,(i+1):end);
                yT = baseLine(plotsubindexT);
                yN = baseLineZ(plotsubindexN);
                ydiffTrue = yT(1:(end-i))-yN((i+1):end);
                improvedX = dataX(plotsubindexT(1:(end-i)),1);
                
                subplot(round(length(estW)/2),round(length(estbbx)/2),(i-1)*round(length(estbbx)/2)+j);
                plot(improvedX,ydiffTrue,'k-');hold on;
                plot(improvedX,nanmean(ydiff,1),'b-');
                scatter(improvedX,quantile(ydiff,0.975,1),'+');
                scatter(improvedX,quantile(ydiff,0.025,1),'o');
                ylim([2*min(min(nanmean(ydiff,1)),min(ydiffTrue)),2*max(max(nanmean(ydiff,1)),max(ydiffTrue))]);
                if (i == 1 && j==1)
                    title(num2str(bands));hold off;
                else
                    title(num2str(bbx));hold off;
                end
                %     title(num2str([bands(1,:);bands(2,:)]));hold off;
            end
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_Improve',num2str(optionsPara(1)),'_',regexprep(num2str([optionsPara(2),optionsPara(3)]),'\s+','_'),'.png')],'png')
        %         close(plotT);
    case 'bias'
        plotT = figure('visible',displayType,'Position', get(0, 'Screensize'));
        bbxallposs = unique(dataX(:,2:end),'rows');
        totalPlotN = length(bbxallposs);
        for i = 1:totalPlotN
            subplot(ceil(sqrt(totalPlotN)),ceil(sqrt(totalPlotN)),i)
            bbx = bbxallposs(i);
            plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows'));
            plot(dataX(plotsubindexT,1),baseLine(plotsubindexT)*0,'-');hold on;
            scatter(dataX(plotsubindexT,1),mean(dataY(:,plotsubindexT),1)'-baseLine(plotsubindexT),'filled');
            scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.975,1)'-baseLine(plotsubindexT),'+');
            scatter(dataX(plotsubindexT,1),quantile(dataY(:,plotsubindexT),0.025,1)'-baseLine(plotsubindexT),'o');
            ylim([-10 max(100,1.2*max(baseLine(plotsubindexT)))]);
            if (i == 1)
                title(num2str(bands));hold off;
            else
                title(num2str(bbx));hold off;
            end
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_N',num2str(optionsPara(1)),'_',regexprep(num2str([optionsPara(2),optionsPara(3)]),'\s+','_'),'.png')],'png')
        %         close(plotT);
    case 'biasChange'
        bbxallposs = unique(dataX(:,2:end),'rows');
        bbxallposs = bbxallposs(1:round(length(bbxallposs(:,1))/9):length(bbxallposs(:,1)),:);
        numOfRowOfPlot = size(dataX,2)-1;
        totalPlotN = length(bbxallposs);
        colors = jet(size(dataY,1));
        plotT = figure('visible',displayType, 'PaperPosition', [0 0 size(dataY,1)*9 9]);
        upper = max(max(quantile(dataY,0.975,1)));
        lower = min(min(quantile(dataY,0.025,1)));
        for ii = 1:size(dataY,1)
%             subplot(numOfRowOfPlot,size(dataY,1),ii);
            subtightplot(numOfRowOfPlot, size(dataY,1), ii, [0.01 0.01], [0.10 0.1], [0.01 0.01]);
            for i = 1:totalPlotN
                bbx = bbxallposs(i,:);
                plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows'));
                plot(dataX(plotsubindexT,1),baseLine(plotsubindexT)*0,'k-');hold on;
                lineMean(i) = plot(dataX(plotsubindexT,1),mean(squeeze(dataY(ii,:,plotsubindexT)),1)-baseLine(plotsubindexT)'...
                    ,'.','MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',30);
                plot(dataX(plotsubindexT,1),quantile(squeeze(dataY(ii,:,plotsubindexT)),0.975,1)-baseLine(plotsubindexT)'...
                    ,'+','MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',10);
                plot(dataX(plotsubindexT,1),quantile(squeeze(dataY(ii,:,plotsubindexT)),0.025,1)-baseLine(plotsubindexT)'...
                    ,'*','MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerSize',10);
            end
            ylim([min(-2,lower) max([10,upper])]);
            xlabel(num2str([optionsPara(2),optionsPara(2+ii)]));
            legend(lineMean,num2str(bbxallposs));
            if (ii == 1)
                title(num2str(bands));hold off;
            else
                title(num2str(optionsPara(2+ii)));hold off;
            end
%         set(gca, 'LooseInset', [0,0,0,0]);
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_N',num2str(optionsPara(1)),'_',num2str(optionsPara(2)),'.png')],'png')
%                 close(plotT);
    case 'closerChange'
        bbxallposs = unique(dataX(:,2:end),'rows');
        numOfRowOfPlot = size(dataX,2)-1;
        totalPlotN = length(bbxallposs);
        colors = jet(totalPlotN);
        plotT = figure('visible',displayType, 'PaperPosition', [0 0 totalPlotN*4 4]);
        for ii = 1:totalPlotN
%             subplot(numOfRowOfPlot,size(dataY,1),ii);
            subtightplot(numOfRowOfPlot, totalPlotN, ii, [0.01 0.01], [0.10 0.1], [0.01 0.01]);
            bbx = bbxallposs(ii);
            plotsubindexT = find(ismember(dataX(:,2:end), bbx,'rows') & dataX(:,1)==optionsPara(3));
            plot(optionsPara(4:end),optionsPara(4:end)*0,'k-');hold on;
            lineMean(ii) = plot(optionsPara(4:end),mean(squeeze(dataY(:,:,plotsubindexT)),2)-baseLine(plotsubindexT),...
                '.','MarkerEdgeColor',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',30);
            plot(optionsPara(4:end),quantile(squeeze(dataY(:,:,plotsubindexT)),0.975,2)-baseLine(plotsubindexT),...
                '+','MarkerEdgeColor',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',10);
            plot(optionsPara(4:end),quantile(squeeze(dataY(:,:,plotsubindexT)),0.025,2)-baseLine(plotsubindexT),...
                '*','MarkerEdgeColor',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',10);
            upper = max(max(quantile(squeeze(dataY(:,:,plotsubindexT)),0.975,2)));
            lower = min(min(quantile(squeeze(dataY(:,:,plotsubindexT)),0.025,2)));
            ylim([min(-2,lower) max([10,upper])]);
            xlabel(num2str([optionsPara(2),optionsPara(3),bbx]));
            legend(lineMean,num2str(bbxallposs));
            if (ii == 1)
                title(num2str(bands));hold off;
            else
                title(num2str(bbx));hold off;
            end
%         set(gca, 'LooseInset', [0,0,0,0]);
        end
        saveas(plotT,[pwd,strcat('\plots\',saveName,'_N',num2str(optionsPara(1)),'_',num2str(optionsPara(2)),'.png')],'png')
%                 close(plotT);
end
end
