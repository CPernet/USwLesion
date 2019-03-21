function [Perf_data, Perf_ranked_data, Perf_adjdata, Perf_ranked_adjdata,cluster_labels] = stats_analysis(csv_folder)

% FORMAT  [Perf_data, Perf_ranked_data, Perf_adjdata, Perf_ranked_adjdata,cluster_labels] = stats_analysis(csv_folder)
%
% INPUT csv_folder is where all the results saved as csv are located
%       e.g. csv_folder = 'C:\Users\s1835343\mri_stuff\spm12\toolbox\USwLesion\validation\BT_analysis_pipe';
%
% OUTPUT Perf_data a 5*48*48 binary matrix indicating which metric is
%                  significantly different from the other
%        Perf_ranked_data a 5*48 matrix of the median rank of each metric
%        --> 5 in dim 1 because Dice was removed after checking the adjusted data
%        --> 48 in dim 2 because thresholding 2 was removed after checking the adjusted data
%        Perf_adjdata a 6*60 matrix indicating improvement or worsening in performence
%        Perf_ranked_adjdata a 6*60 matrix indicating if ranking improved or worsend 
%        cluster_labels a cell array of labels for the dendrograms thresholded at 
%                       [2 4 3 6 8 12 24] corresponding to our manipulations

% 1 - checking if there are significant increase or decrease in perf 
% 2 - check how data cluster

for m=2:-1:1 % backward to start with 60 param on difference then 48 of raw
    for d = 1:6
        if d==1
            name = 'mask_volumes.csv'; nname = 'Volumes';
        elseif d==2
            name = 'mJ_results.csv'; nname = 'Jaccard';
        elseif d==3
            name = 'Dice_results.csv'; nname = 'Dice';
        elseif d==4
            name = 'mHd_results.csv'; nname = 'Hausdorff';
        elseif d== 5
            name = 'overlap_kappa_results.csv'; nname = 'Kappa';
        elseif d==6
            name = 'overlap_mcc_results.csv'; nname = 'Matthew Corr';
        end
        
        IMP   = importdata([csv_folder filesep name]);
        data  = IMP.data;
        label = IMP.colheaders;
        clear IMP
        
        if m == 1
            data(:,13:24) = []; % remove threshold 2 because we know from looking at the difference it performs badly
            label(13:24)  = [];
        elseif m==2
            if d == 1
                IMP  = importdata([csv_folder filesep 'ground_truth_volumes.csv']);
                data = data - repmat(IMP.data,1,60);
            else
                IMP = importdata([csv_folder filesep 'baseline_measures.csv']);
                if d==2 % Jaccard
                    data = data - repmat([repmat(IMP.data(:,1),[1 6]) repmat(IMP.data(:,2),[1 6]) ],[1,5]);
                elseif d==3 % Dice
                    data = data - repmat([repmat(IMP.data(:,5),[1 6]) repmat(IMP.data(:,6),[1 6]) ],[1,5]);
                elseif d==4 % Hausdorff
                    data = repmat([repmat(IMP.data(:,3),[1 6]) repmat(IMP.data(:,4),[1 6]) ],[1,5]) - data;
                elseif d==5 % Kappa
                    data = data - repmat([repmat(IMP.data(:,9),[1 6]) repmat(IMP.data(:,10),[1 6]) ],[1,5]);
                elseif d== 6 % Matthew Corr
                    data = data - repmat([repmat(IMP.data(:,7),[1 6]) repmat(IMP.data(:,8),[1 6]) ],[1,5]);
                end
                clear IMP
            end
        end
        
        % data viz and ranking
        % ---------------------
        figure
        subplot(2,1,1);
        [~,HDI] = rst_data_plot(data,'estimator','trimmed mean','newfig','no');
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
        
        med = sort(med);
        prob_coverage = 95/100;
        upper_centile = floor(prob_coverage*size(med,1)); % upper bound
        nCIs = size(med,1) - upper_centile;
        ci = 1:nCIs; ciWidth = med(ci+upper_centile) - med(ci); % all centile distances
        [~,index]=find(ciWidth == min(ciWidth)); % densest centile
        if length(index) > 1; index = index(1); end % many similar values
        RHDI(d,1,:) = med(index,:);
        RHDI(d,2,:) = med(index+upper_centile,:);
        rst_boxplot(med); hold on
        plot([1:length(RHDI)],squeeze(RHDI(d,1,:)),'--r','LineWidth',2);
        plot([1:length(RHDI)],squeeze(RHDI(d,2,:)),'--r','LineWidth',2);
        if m == 1
            title([nname ' Ranked Raw Data']) 
        else
            title([nname ' Ranked Adjusted Data'])
        end
        clear med index nCI ci ciWidth
        
        % check if HDI/RHDI overlaps for statistical difference
        if m == 1
            if d == 1 % create array to fill
                Perf_data        = NaN(6,48,48);
                Perf_ranked_data = zeros(6,48);
            end
            
            % check which params are higher than others
            for p=1:48
                Perf_data(d,p,:) = HDI(1,:)>HDI(2,p);
            end
            
            [~,index] = sort(data','descend');
            Perf_ranked_data(d,:) = rst_hd(index',0.5);
            
            if d == 6 % when all is analyzed, make a summary figure 
               Perf_data(3,:,:)        = []; % remove Dice as we know from the difference analysis it performes badly
               Perf_ranked_data(3,:,:) = []; 
               RHDI(3,:,:)             = [];
               
               figure; subplot(1,2,1); imagesc(squeeze(sum(Perf_data,1)));
               title('frequency of HDI differring')
               subplot(1,2,2); plot([1:48],Perf_ranked_data); grid on
               hold on; plot([1:48], mean(Perf_ranked_data),'k','LineWidth',2);
               plot([1:48], squeeze(mean(RHDI)),'-k','LineWidth',2);
               plot([1:48], squeeze(mean(RHDI)),'-k','LineWidth',2);
               title('Median rank per metric (and mean across metrics)'); 
               clear RHDI HDI
            end
            
            % compute similarity and clustering
            dist = pdist(data','euclidean');
            figure; subplot(1,2,1); 
            imagesc(squareform(dist)); title(nname)
            L = linkage(dist,'average');
            subplot(1,2,2);
            dendrogram(L,'Orientation','left');
            C = cophenet(L,dist); title(sprintf('clustering coef %g',C))
            
            % check clusters at different levels (eg do VOI1 VOI2 get separated?)
            % 2 roi 4 thresholdings 3 gaussians 2 tissues
            threshold = [2 4 3 6 8 12 24];
            for t=1:length(threshold)
                T(:,t) = cluster(L,'maxclust',threshold(t));
                for l=1:length(unique(T(:,t)))
                    cluster_labels{d}{t,l} = (label(T(:,t)==l))';
                end
            end
            
        end
        
        % check if algorithm has improved or worsened similarity of mask to ground truth
        % improved if lower bound > 0; worsened if upper bound < 0
        if m == 2
            if d == 1 % create array to fill
                Perf_adjdata        = zeros(5,60);
                Perf_ranked_adjdata = zeros(5,60);
            end
            Perf_adjdata(d,HDI(1,:)>0)           = 1;
            Perf_adjdata(d,HDI(2,:)<0)           = -1;
            Perf_ranked_adjdata(d,RHDI(d,1,:)>0) = 1;
            Perf_ranked_adjdata(d,RHDI(d,2,:)<0) = -1;
            
            if d == 6 
                figure; subplot(2,1,1); imagesc(Perf_adjdata);
                subplot(2,1,2); imagesc(Perf_ranked_adjdata);
                disp('adjusted data')
                disp(['performing better ' num2str(find(sum(Perf_adjdata == 1)))])
                disp(['performing worst ' num2str(find(sum(Perf_adjdata == -1)))])
                disp(['showing both better and worst perf ' num2str(intersect(find(sum(Perf_adjdata == 1)),find(sum(Perf_adjdata == -1))))])
                clear RHDI HDI
            end
        end
        
    end
end
