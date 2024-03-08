%            Date：      20230102
%            Author:     LvQIuyv
% --------introduction-------------------------------------------------------------------------------
%           This program implement predictive models based on brain
%           Input :     1:  all subject functonal connectivity matrix  csv files  at the 20th line of code :folder = "E:/LEMON/CPM_shen/CorrMatrix_40";
%           
% 
%                          2:   behavioral data:the 35rd line of code and
%           the 36th line of code are respectively the lines where all the behavioral data csv files and the behavioral data to be predicted are located
% 
%           Output:      1. p-value and r-value of predicted effect and  csv file of predicted values of all behavioral data
%                             2. 
% 
% 
% 
%%------------------------------------------------------------------------------------------------------
clear;
clc;
% 指定路径（唯一需要需要修改的地方）
folder  = "G:\Graduation Project\OneDrive - hunnu\STUDY2\Matrix\Matrix_FishZ\TG_268_F";
behavdatapath = "G:/Graduation Project/OneDrive - hunnu/STUDY2/Behav/WM.csv";
columnposition = 4;
outputOptiR =  'G:/Graduation Project/OneDrive - hunnu/STUDY2/OC/SWMS/TG_268/1_Optimal_r_value_LOOCV.csv';

% -------- readMatrix/behav -----------

files = dir(fullfile(folder, '*.csv'));
% Initialize the 3D matrix
data = [ ];

% Loop over all csv files
for i = 1:length(files)
    % Load each csv file into a matrix
    file = fullfile(folder, files(i).name);
    csv_data = csvread(file);
    
    % Add the matrix to the 3D matrix along the third dimension
    data(:,:,i) = csv_data;
end

behavdata = readtable(behavdatapath);
BAS_Drive = table2array(behavdata(:,columnposition));

% ------------ INPUTS -------------------

all_mats  = data;
all_behav = BAS_Drive;



no_sub = size(all_mats,3);
no_node = size(all_mats,1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred = zeros(no_sub,1);

% threshvalues = 0.0001 : 0.0001:0.0005;
% threshold for feature selection
% thresh = 0.01;
Optimal_r = [];
% ---------------------------------------
for thresh = 0.0001 : 0.0001:0.05
fprintf('\n 当前阈值为 # %10.5f\n',thresh);

    for leftout = 1:no_sub
    fprintf('\n Leaving out subj # %6.3f\n',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    
    edge_no = size(train_vcts,1);
    r_mat = zeros(1, edge_no);
    p_mat = zeros(1, edge_no);
%     
%     Use different methods to correlate all edges with behavior:
% 
%     
%           Method1:    correlate all edges with behavior using robust regression
% 
%           for edge_i = 1: edge_no;
%           [~, stats] = robustfit(train_vcts(edge_i,:)', train_behav);
%           cur_t = stats.t(2);
%           r_mat(edge_i) = sign(cur_t)*sqrt(cur_t^2/(no_sub-1-2+cur_t^2));
%           p_mat(edge_i) = 2*(1-tcdf(abs(cur_t), no_sub-1-2));  %two tailed
%           end
% 
% 
%           Method2:    correlate all edges with behavior using partial correlation
%           [r_mat, p_mat] = partialcorr(train_vcts', train_behav, age);
%     
%        
%           Method3:    correlate all edges with behavior using rank correlation
             [r_mat, p_mat] = corr(train_vcts', train_behav, 'type', 'Pearson');
    
        
            r_mat = reshape(r_mat,no_node,no_node);
            p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks 
    pos_mask = zeros(no_node, no_node);
    neg_mask = zeros(no_node, no_node);
    
    
    pos_edge = find( r_mat >0 & p_mat < thresh);
    neg_edge = find( r_mat <0 & p_mat < thresh);
    
    pos_mask(pos_edge) = 1;
    neg_mask(neg_edge) = 1;
    
    
%     %-----------------sigmoidal weighting---------------------------%
%     pos_edges = find(r_mat > 0 );
%     neg_edges = find(r_mat < 0 );
%     
%     % covert p threshold to r threshold
%     T = tinv(thresh/2, no_sub-1-2);
%     R = sqrt(T^2/(no_sub-1-2+T^2));
%     
%     % create a weighted mask using sigmoidal function
%     % weight = 0.5, when correlation = R/3;
%     % weight = 0.88, when correlation = R;
%     pos_mask(pos_edges) = sigmf( r_mat(pos_edges), [3/R, R/3]);
%     neg_mask(neg_edges) = sigmf( r_mat(neg_edges), [-3/R, R/3]);
%     %---------------sigmoidal weighting-----------------------------%
    
    % get sum of all edges in TRAIN subs (divide by 2 to control for the
    % fact that matrices are symmetric)
    
        train_sumpos = zeros(no_sub-1,1);
        train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos);
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end
    
    % build model on TRAIN subs
    % combining both postive and negative features
    b = regress(train_behav, [train_sumpos, train_sumneg, ones(no_sub-1,1)]);
    
    %%分别拟合pos和neg网络
    fit_pos = polyfit(train_sumpos,train_behav,1);
    fit_neg = polyfit(train_sumneg,train_behav,1);
    % run model on TEST sub
    
    test_mat = all_mats(:,:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
    
    behav_pred(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    
    end


        % compare predicted and observed scores

        [R_pos, P_pos] = corr(behav_pred_pos,all_behav);
        [R_neg, P_neg] = corr(behav_pred_neg,all_behav);

        [R, p] = corr(behav_pred, all_behav);
        format shortg
        R = round(R, 5);
        p = round(p, 5);
        format

        output = [thresh,R_pos, P_pos,R_neg, P_neg,R,p];
        numRows = size(Optimal_r, 1);
        index = numRows + 1;
        Optimal_r(index, :) = output;

end
%%%% 1  find the max r_value
% % % %  2  Find the corresponding thresh
MaxVals = max(max(Optimal_r(:,[2 4 6])));
[row,col] = find(Optimal_r == MaxVals);
MaxThresh = Optimal_r(row,1)

fprintf('When the threshold is %0.4f, the model is optimal, where the rho at this time is %0.4f.\n', MaxThresh, MaxVals);

columnname = {'thresh','R_pos', 'P_pos','R_neg', 'P_neg','R','p'};
threshdata = table(Optimal_r(:,1),Optimal_r(:,2),Optimal_r(:,3),Optimal_r(:,4),Optimal_r(:,5),Optimal_r(:,6),Optimal_r(:,7),'VariableName',columnname);
writetable(threshdata,outputOptiR);
%-------------------------------------------------------------------------------------------------------------------------
% Write a matrix containing the r-values and p-values of the entire network, positive network, and negative network:
% 
%           two columns:   r_value、p_value     
%           three rows:      positive net,     negative net,     all net 
%           
% ----------------------------------------------------------------------------------------
% Net_Name  = {'Pos_Net';'Neg_Net';'Whole_Net'};
% p_value = [P_pos;P_neg;p];
% r_value = [R_pos;R_neg;R];
% column = {'Net_Name','r_value','p_value'};
% RPdata = table(Net_Name,r_value,p_value,'VariableNames',column);
% writetable(RPdata,'./Result/1_RP_data.csv');
% 


