function [Perf_data, Perf_ranked_data, Perf_adjdata, Perf_ranked_adjdata] = stats_analysis(csv_folder)

% csv_folder = '/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/toolbox/USwLesion/validation/BT_analysis_pipe';

for m=1:2
    for d = 1:5
        if d==1
            name = ['Dice_results.csv']; nname = 'Dice';
        elseif d==2
            name = ['mHd_results.csv']; nname = 'Hausdorff';
        elseif d==3
            name = ['mJ_results.csv']; nname = 'Jaccard';
        elseif d==4
            name = ['overlap_kappa_results.csv']; nname = 'Kappa';
        elseif d== 5
            name = ['overlap_mcc_results.csv']; nname = 'Matthew Corr';
        end
        
        IMP = importdata([csv_folder filesep name]);
        data = IMP.data;
        label = IMP.colheaders;
        clear IMP
        
        if m==2
            IMP = importdata([csv_folder filesep 'baseline_measures.csv']);
            if d==1 % Dice
                data = data - repmat([repmat(IMP.data(:,5),[1 6]) repmat(IMP.data(:,6),[1 6]) ],[1,5]);
            elseif d==2 % Hausdorff
                data = data - repmat([repmat(IMP.data(:,3),[1 6]) repmat(IMP.data(:,4),[1 6]) ],[1,5]);
            elseif d==3 % Jaccard
                data = data - repmat([repmat(IMP.data(:,1),[1 6]) repmat(IMP.data(:,2),[1 6]) ],[1,5]);
            elseif d==4 % Kappa
                data = data - repmat([repmat(IMP.data(:,9),[1 6]) repmat(IMP.data(:,10),[1 6]) ],[1,5]);
            elseif d== 5 % Matthew Corr
                data = data - repmat([repmat(IMP.data(:,7),[1 6]) repmat(IMP.data(:,8),[1 6]) ],[1,5]);
            end
            clear IMP
        end
        
        % data viz and ranking
        
        figure
        subplot(2,1,1); 
        [est,HDI] = rst_data_plot(data,'estimator','trimmed mean','newfig','no');
        if m == 1
            title([nname ' Raw Data']);
        else
            title([nname ' Adjusted Data']);
        end
        % bootstrap the data and rank
        subplot(2,1,2);
        [o,n] = size(data);
        med = NaN(599,n);
        for boot=1:599 % bootstrap loop
            theta = exprnd(1,[n,1]);
            weigths = theta ./ repmat(sum(theta),n,1);
            resample = (datasample(data',n,'Replace',true,'Weights',weigths));
            [~,index] = sort(resample,'descend');
            med(boot,:) = rst_hd(index',0.5);
        end
        rst_boxplot(med); 
        if m == 1
            title([nname ' Ranked Data']) 
        else
            title([nname ' Adjusted Ranked Data'])
        end
        med = sort(med);
        prob_coverage = 95/100;
        upper_centile = floor(prob_coverage*size(med,1)); % upper bound
        nCIs = size(med,1) - upper_centile;
        ci = 1:nCIs; ciWidth = med(ci+upper_centile) - med(ci); % all centile distances
        [~,index]=find(ciWidth == min(ciWidth)); % densest centile
        if length(index) > 1; index = index(1); end % many similar values
        RHDI(1,:) = med(index,:);
        RHDI(2,:) = med(index+upper_centile,:);
        
        % check if HDI/RHDI overlaps for statistical difference 
        if m == 1
            if d == 1 % create array to fill
                Perf_data        = NaN(5,60,60);
                Perf_ranked_data = zeros(5,60);
            end
            
            % check which params are higher than others
            for p=1:60
                Perf_data(d,p,:) = HDI(1,:)>HDI(2,p);
            end
            
            [~,index] = sort(data','descend');
            Perf_ranked_data(d,:) = rst_hd(index',0.5);

            if d == 5 
               figure; subplot(1,2,1); imagesc(squeeze(sum(Perf_data,1)));
               title('frequency of HDI differring')
               subplot(1,2,2); plot([1:60],Perf_ranked_data); grid on
               hold on; plot([1:60], mean(Perf_ranked_data),'-k','LineWidth',2);
               title('Median rank per metric (and mean across metrics)'); 
            end
            
            
            % compute similarity and clustering
            dist = pdist(data','euclidean');
            figure; subplot(1,2,1); 
            imagesc(squareform(dist)); title(nname)
            L = linkage(dist,'average');
            subplot(1,2,2);
            [~,T] = dendrogram(L,'Orientation','left');
            C = cophenet(L,dist); title(sprintf('clustering coef %g',C))
            
            % check clusters at different levels (eg do VOI1 VOI2 get separated?)
            
            % T = cluster(L,'cutoff',?);
            
        end
        
        % check if algorithm has improved or worsened similarity of mask to ground truth
        % improved if lower bound > 0; worsened if upper bound < 0
        if m == 2
            if d == 1 % create array to fill
                Perf_adjdata        = zeros(5,60);
                Perf_ranked_adjdata = zeros(5,60);
            end
            Perf_adjdata(d,HDI(1,:)>0)        = 1;
            Perf_adjdata(d,HDI(2,:)<0)        = -1;
            Perf_ranked_adjdata(d,RHDI(1,:)>0) = 1;
            Perf_ranked_adjdata(d,RHDI(2,:)<0) = -1;
            
            if d == 5 
                figure; subplot(2,1,1); imagesc(Perf_adjdata);
                subplot(2,1,2); imagesc(Perf_ranked_adjdata);
                disp('adjusted data')
                disp(['performing better ' num2str(find(sum(Perf_adjdata == 1)))])
                disp(['performing worst ' num2str(find(sum(Perf_adjdata == -1)))])
                disp(['showing both better and worst perf ' num2str(intersect(find(sum(Perf_adjdata == 1)),find(sum(Perf_adjdata == -1))))])
            end
        end
        
    end
end
