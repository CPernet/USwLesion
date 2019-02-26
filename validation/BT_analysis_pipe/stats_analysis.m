function stats_analysis(csv_folder);

% csv_folder = '/home/cpernet/Documents/MATLAB/MyToolboxes/SPM/spm12/toolbox/USwLesion/validation/BT_analysis_pipe';


for d = 1:5
    if d==1
        name = ['Dice_results.csv']; nname = 'Dice';
    elseif d==2
        name = ['mHd_results.csv']; nname = 'Hansdorff';
    elseif d==3
        name = ['mj_results.csv']; nname = 'Jaccard';
    elseif d==4
           name = ['overlap_kappa_results.csv']; nname = 'Kappa';
    elseif d== 5
        name = ['overalap_mcc_results.csv']; nname = 'Matthew Corr';
    end
    
    IMP = importdata([csv_folder filesep name]);
    data = IMP.data;
    label = IMP.colheaders;
    clear IMP
    
% ranking the data in ascending order of similarity for each patient. for each overlap metric, sd is the sorted data and i is the index.

   figure
   subplot(2,1,1); [est,HDI] = rst_data_plot(data,'estimator','trimmed mean','newfig','no'); 
   title([nname ' Raw Data']);  
   % boostrap the data and rank
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
   rst_boxplot(med); title([nname ' Ranked Data'])
   med = sort(med);
   prob_coverage = 95/100;
   upper_centile = floor(prob_coverage*size(med,1)); % upper bound
   nCIs = size(med,1) - upper_centile;
   ci = 1:nCIs; ciWidth = med(ci+upper_centile) - med(ci); % all centile distances
   [~,index]=find(ciWidth == min(ciWidth)); % densest centile
   if length(index) > 1; index = index(1); end % many similar values
   RHDI(1,:) = med(index,:);
   RHDI(2,:) = med(index+upper_centile,:);
   
   
% clustering the data and visualising in a dendrogram 

dist = pdist(data','euclidean');
figure; subplot(1,2,1); imagesc(squareform(dist))
L = linkage(dist,'average');
subplot(1,2,2); 
[H,T,cluster_labels] = dendrogram(L,'Orientation','left','Labels' ,label);
C = cophenet(L,dist);


end

