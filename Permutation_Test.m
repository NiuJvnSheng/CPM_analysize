function Permutation_Test(PATH_Matrix, PATH_Behav, ncol, thresh_f, trial, output)

    % 1 导入matrix数据
    folder = PATH_Matrix;
    files = dir(fullfile(folder, '*.csv'));

    % Initialize the 3D matrix
    data = [];
    for i = 1:length(files)
        % Load each csv file into a matrix
        file = fullfile(folder, files(i).name);
        csv_data = csvread(file);

        % Add the matrix to the 3D matrix along the third dimension
        data(:,:,i) = csv_data;
    end
    % 2 导入行为学数据
    behavdata = readtable(PATH_Behav);
    BAS_Drive = table2array(behavdata(:,ncol));
    all_mats  = data;
    all_behav = BAS_Drive;
    % 3.1 定义被试人数和置换检验的次数
    no_sub = size(all_mats,3);
    % 3.2 定义置换检验的次数
    times = trial;
    % 3.2 定义置换检验的for list
    no_iterations   = times+1;

    [true_prediction_r_pos, true_prediction_r_neg,true_prediction_r_all] = predict_behavior(all_mats, all_behav,thresh_f);
    prediction_r    = zeros(no_iterations,3);
    prediction_r(1,1) = true_prediction_r_pos;
    prediction_r(1,2) = true_prediction_r_neg;
    prediction_r(1,3) = true_prediction_r_all;

    parfor it=2:no_iterations
        fprintf('\n Performing iteration %d out of %d\n', it, no_iterations)
        new_behav = all_behav(randperm(no_sub));
        [temp_pred_r_pos, temp_pred_r_neg, temp_pred_r_all] = predict_behavior(all_mats, new_behav,thresh_f);
        prediction_r(it,:) = [temp_pred_r_pos, temp_pred_r_neg, temp_pred_r_all];
    end

    Pos_R  = prediction_r(:,1);
    Neg_R = prediction_r(:,2);
    All_R = prediction_r(:,3);
    column = {'Pos_R','Neg_R','All_R'};
    Permutation_Data = table(Pos_R,Neg_R,All_R,'VariableNames',column);
    writetable(Permutation_Data, output);

end