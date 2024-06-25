%% Init
close all;
clear; clc;

% Paths
%dataCNDFolder = '.\outputs\64Hz\CND\';
%dataStimFolder = '.\outputs\64Hz\features\';
% Changing to Mac path
%dataCNDFolder = './eegProject/datasets/SLdata/1-30Hz/Mastoids/';
%dataStimFolder = './outputs/30Hz/features/';

% June 19, 2024 new output path
dataCNDFolder = './outputs_new/64Hz/CND/';
dataStimFolder = './outputs_new/64Hz/features/';

% Add other directories to path
% addpath ..\..\MATLAB\cnsp_utils
% addpath ..\..\MATLAB\cnsp_utils\cnd
addpath ./CNSP-resources/CNSP/libs/cnsp_utils/
addpath ./CNSP-resources/CNSP/libs/cnsp_utils/cnd/

%% General parameters
% Number of trials for each condition
NTRIALS = 14;
% Downsampled frequency
fs_down = 30;
% Conditions: video and reverse
conditions = {'V', 'R'}; 
feature_names = {'envelope1', 'envelope2'};
additionalDetails = {'envelope1', 'envelope2'};
% Features that need to be transposed (flipped) - because of those reveresed ones
feats2transpose = {'envelope1', 'envelope2'};

%% Generate data stim for each condition
for condition = conditions
    cond = condition{1};
    
    stim = struct();
    stim.fs = fs_down;
    stim.condition = cond;
    for idx_feat = 1:length(feature_names)
        feature_name = feature_names{idx_feat};
%         stim_file_names = dir([dataStimFolder,feature_name,'\',cond,'*','.mat']);

        stim.names{idx_feat} = feature_name;
        stim.additionalDetails{idx_feat} = additionalDetails{idx_feat};
        for idx_stim = 1:NTRIALS
            if idx_stim > 9
                stim_name = [cond, num2str(idx_stim), '.mat'];
            else
                stim_name = [cond, '0', num2str(idx_stim), '.mat'];
            end
            stim_path = [dataStimFolder,feature_name,'/',stim_name];
            disp(stim_path)
            stim.audioFiles{idx_stim} = stim_name;

            if ~exist(stim_path, 'file')
                stim.data{idx_feat, idx_stim} = [];
                disp([stim_name, ' does not exist'])
            else
                % Save feature data in cell array
                load(stim_path,'feature')
                feature(isnan(feature)) = 0;
                if ~any(strcmp(feats2transpose,feature_name))
                    stim.data{idx_feat, idx_stim} = feature;
                else
                    stim.data{idx_feat, idx_stim} = feature';
                end
            end

        end

    end

    outputFilename = [dataCNDFolder, cond, '/dataStim.mat'];
    disp(['Saving data stim for condition ', cond])
    save(outputFilename,'stim')
end
 
