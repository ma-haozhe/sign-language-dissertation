%% Init
close all;
clear; clc;

%% Parameters
% Select EEG data to use
datafolder = '.\outputs\64Hz\CND\';
TRFfolder = '.\outputs\64Hz\';
nSubs = 2;
conditions = {'V', 'R'}; 
reRefType = 'Mastoids'; 
pre = '1-30Hz'; 

sfreq = 64;
remove_start_and_end = 1;
rm_seconds = 0.5; % s
rm_samples = rm_seconds*sfreq;

% Select features to use
dataStim = '\dataStim.mat';
feature_name = 'A2'; 
feature_idxs = [2]; % IdyomTony=21,Phentr=36,Pht=35
shuffle_idxs = [1];  
onsets_idxs  = [1]; % if shuffle_type='shuffle_values' %noteonsets=29,phonsets=32,f0=11
shuffle_mode = 1;
shuffle_type = 'shuffle';  %'shuffle'; % 'flip'; %'shuffle_values';
idx_feat_int = 1;
xlim_values = [0, 400];
margin = 50;
t1 = 30;
t2 = 94;
t3 = 500;

% Model hyperparameters
chan = 38;   % Fz (only for decoding results display)
Dir = 1;     % forward=1, backward=-1 
tmin = xlim_values(1) - margin;
tmax = xlim_values(2) + margin;
lamda_idx = -4:1:6;
lambda_vals = 10.^lamda_idx; % 1000
nlambda = numel(lambda_vals);

% Output path
tauminstr = num2str(tmin);
taumaxstr = num2str(tmax);

% Output path
if Dir == 1
    fold = 'TRFs';
else
    fold = 'Decoders';
end

output_path = [TRFfolder, fold,'/', pre,'/', tauminstr, '_', taumaxstr, '/', reRefType, '/'];
mkdir(output_path);
addpath(output_path);

for condition = conditions
    cond = condition{1};

    %% Load features
    disp('Loading data...')
    stimPathname = [datafolder, cond, dataStim];
    disp(['Loading stimulus data: ', stimPathname])
    load(stimPathname,'stim');

    %% Nasted crossvalidation
    for sub = 1:nSubs
        %% Load preprocessed EEG
        eegPreFilename = [datafolder,cond,'/',pre,'/',reRefType,'\pre_dataSub', num2str(sub),'.mat'];
        load(eegPreFilename,'eeg')
        FS_EEG = eeg.fs;

        %% Deal with empty trials and do shuffeling
        data = eeg.data;
        i = 0;
        for tr = 1:length(data)
            tmpEeg = data{tr};
            if ~isempty(tmpEeg)
                i = i + 1;
                response{i} = tmpEeg; % clear tmpEeg;
                
                tmpEnv = [];
                tmpShu = [];
                for j = feature_idxs    
                    if shuffle_idxs(feature_idxs==j)
                        shuffled_feature = stim.data{j, tr};
                        if strcmp(shuffle_type, 'flip')
                            shuffled_feature = flip(shuffled_feature);
                        elseif strcmp(shuffle_type, 'shuffle')
                            shuffled_feature = shuffle(shuffled_feature);
                        elseif strcmp(shuffle_type, 'shuffle_values')
                            idx_feat_onset = onsets_idxs(feature_idxs==j);
                            onsets = stim.data{idx_feat_onset, tr};
                            onsets_loc = onsets > 0;
                            values = shuffled_feature(onsets_loc);
                            values = shuffle(values);
                            shuffled_feature2 = zeros(size(shuffled_feature));
                            shuffled_feature2(onsets_loc) = values;
                            shuffled_feature = shuffled_feature2;
                        end
                        tmpShu = [tmpShu, shuffled_feature];
                    else
                        tmpShu = [tmpEnv, stim.data{j, tr}];
                    end
                    tmpEnv = [tmpEnv, stim.data{j, tr}];
                end
                features{i} = tmpEnv; clear tmpEnv; 
                shuffele{i} = tmpShu; clear tmpShu;
            end
        end
        n_training_trial = i;

        %% Making sure that stim and neural data have the same length
        for tr = 1:length(features)
            envLen = size(features{tr},1);
            eegLen = size(response{tr},1);
            minLen = min(envLen,eegLen);
            features{tr} = double(features{tr}(1:minLen,:));
            shuffele{tr} = double(shuffele{tr}(1:minLen,:));
            response{tr} = double(response{tr}(1:minLen,:));
        end

        %% Remove start and end of trial
        if remove_start_and_end
            response = cellfun(@(x) x(rm_samples:size(x, 1)-rm_samples, :),response,'UniformOutput',false);
            features = cellfun(@(x) x(rm_samples:size(x, 1)-rm_samples, :),features,'UniformOutput',false);
            shuffele = cellfun(@(x) x(rm_samples:size(x, 1)-rm_samples, :),shuffele,'UniformOutput',false); 
        end

        %% Normalising EEG data
        eeg_data_mat = cell2mat(response');
        eeg_mean = mean(eeg_data_mat(:));
        eeg_std = std(eeg_data_mat(:));
        response = cellfun(@(x) (x-eeg_mean)/eeg_std,response,'UniformOutput',false);
        
        %% Crossvalidation for subject j
        % Define training and test sets
        stim_train = features; 
        eeg_train = response;
    
        % Run fast cross-validation
        disp('Running cross-validation...')
        cv = mTRFcrossval(stim_train,eeg_train,FS_EEG,Dir,tmin,tmax,lambda_vals);
        
        % Get optimal hyperparameters
        [rmax,idx] = max(mean(mean(cv.r, 1), 3));
        lambda = lambda_vals(idx);
    
        % Train model
        disp('Training model...')
        Smodel = mTRFtrain(stim_train,eeg_train,FS_EEG,Dir,tmin,tmax,lambda,'zeropad',0);
        Smodel.lambda = lambda;

        if Dir == 1
             % Plot CV accuracy and error
            figure()
            sgtitle([cond, ' sub ', num2str(sub), ' nr trials ', num2str(n_training_trial), ' best lambda ', num2str(lambda)])
            subplot(1,2,1)
            errorbar(1:nlambda,mean(cv.r(:,:,chan)),std(cv.r(:,:,chan))/sqrt(numel(stim_train)),'linewidth',2)
            set(gca,'xtick',1:nlambda,'xticklabel',lamda_idx), xlim([0,nlambda+1])
            title(['Best lambda ', num2str(lambda)])
            xlabel('Regularization (1\times10^\lambda)')
            ylabel('Correlation')
            axis square, grid on
            subplot(1,2,2)
            errorbar(1:nlambda,mean(cv.err(:,:,chan)),std(cv.err(:,:,chan))/sqrt(numel(stim_train)),'linewidth',2)
            set(gca,'xtick',1:nlambda,'xticklabel',lamda_idx), xlim([0,nlambda+1])
            title('CV Error')
            xlabel('Regularization (1\times10^\lambda)')
            ylabel('MSE')
            axis square, grid on
    
            % Plot TRF, GFP and TRF weights
            [topo_minValue,idx_min] = min(abs(Smodel.t-tmin));   
            [topo_maxValue,idx_max] = min(abs(Smodel.t-tmax));
            [topo1_val, topo1_idx] = min(abs(Smodel.t - t1));
            [topo2_val, topo2_idx] = min(abs(Smodel.t - t2));
            lim = max(max(abs(Smodel.w(idx_feat_int,idx_min:idx_max,:)),[],3),[],2);
    
            figure()
            sgtitle([cond, ' sub ', num2str(sub), ' nr trials ', num2str(n_training_trial), ' best lambda ', num2str(lambda)])
            subplot(2,2,1)
            plot(Smodel.t,squeeze(Smodel.w(idx_feat_int,:,:)))
            hold on
            yline(0, '-', 'Alpha', 0.9)
            xline(0, '-', 'Alpha', 0.9) 
            xline(t1, '--', 'Alpha', 0.9) 
            xline(t2, '--', 'Alpha', 0.9) 
            xline(t3, '--', 'Alpha', 0.9) 
            xlim(xlim_values)
            title('Temporal Response Function (TRF)')
            xlabel('Time lag (ms)')
            ylabel('Amplitude (a.u.)')
            subplot(2,2,2)
            plot(Smodel.t,std(Smodel.w(idx_feat_int,:,:),[],3))
            hold on
            yline(0, '-', 'Alpha', 0.9)
            xline(0, '-', 'Alpha', 0.9) 
            xline(t1, '--', 'Alpha', 0.9) 
            xline(t2, '--', 'Alpha', 0.9) 
            xline(t3, '--', 'Alpha', 0.9) 
            xlim(xlim_values)
            title('Global Field Power (GFP)')
            xlabel('Time lag (ms)')
            subplot(2,2,3)
            topoplot(Smodel.w(idx_feat_int,topo1_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
            title([num2str(Smodel.t(topo1_idx)),' ms'])
            subplot(2,2,4)
            topoplot(Smodel.w(idx_feat_int,topo2_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
            title([num2str(Smodel.t(topo2_idx)),' ms'])
        
    
    %          % Plot CV accuracy and error
    %         sgtitle([cond, ' subject ', num2str(sub)])
    %         subplot(3,nSubs,1+(sub-1))
    %         errorbar(1:nlambda,mean(cv.r(:,:,chan)),std(cv.r(:,:,chan))/sqrt(numel(stim_train)),'linewidth',2)
    %         set(gca,'xtick',1:nlambda,'xticklabel',lamda_idx), xlim([0,nlambda+1])
    %         title(['Best lambda ', num2str(lambda)])
    %         xlabel('Regularization (1\times10^\lambda)')
    %         ylabel('Correlation')
    %         subplot(3,nSubs,3+(sub-1))
    %         plot(Smodel.t,squeeze(Smodel.w))
    %         hold on
    %         yline(0, '-', 'Alpha', 0.5)
    %         xline(0, '-', 'Alpha', 0.5) 
    %         xline(t1, '--', 'Alpha', 0.9) 
    %         xline(t2, '--', 'Alpha', 0.9) 
    %         xline(t3, '--', 'Alpha', 0.9) 
    %         xlim(xlim_values)
    %         title('Temporal Response Function (TRF)')
    %         xlabel('Time lag (ms)')
    %         ylabel('Amplitude (a.u.)')
    %         subplot(3,nSubs,8+(sub-1))
    %         plot(Smodel.t,std(Smodel.w,[],3))
    %         xlim(xlim_values)
    %         title('Global Field Power (GFP)')
    %         xlabel('Time lag (ms)')
        end

        % Save model and predictions to avarage later
        modelAll(sub) = Smodel; clear Smodel;
        rpredAll(sub, :) = squeeze(mean(cv.r(:, idx, :), 1));
        rpredAll_trials{sub} = squeeze(cv.r(:, idx, :)); clear cv;

        if shuffle_mode
            cv_shu = mTRFcrossval(shuffele,eeg_train,FS_EEG,Dir,tmin,tmax,lambda_vals);
            [rmax_shu,idx_shu] = max(mean(mean(cv_shu.r, 1), 3));
            rpredShu(sub, :) = squeeze(mean(cv_shu.r(:, idx_shu, :), 1)); 
            rpredShu_trials{sub} = squeeze(cv_shu.r(:, idx_shu, :)); clear cv_shu;
        end               
    end

    prename = [feature_name, '_', cond, '_'];
    save([output_path, prename, 'modelAll_cvonly.mat'],'modelAll'); clear modelAll;
    save([output_path, prename, 'rpredAll_cvonly.mat'],'rpredAll'); clear rpredAll;  
    save([output_path, prename, 'rpredAll_trials_cvonly.mat'],'rpredAll_trials'); clear rpredAll_trials; 
    if shuffle_mode
        save([output_path, prename, 'rpredShu_cvonly.mat'],'rpredShu'); clear rpredShu;
        save([output_path, prename, 'rpredShu_trials_cvonly.mat'],'rpredShu_trials'); clear rpredShu_trials;
    end
end

