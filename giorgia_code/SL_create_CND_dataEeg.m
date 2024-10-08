%% Generate CND files from data segmented and preprocessed in Python
%% Do bad channel rejection here as opposed to Python
%% Init
close all;
clear; clc;

% Add other directories to path

% addpath ..\..\MATLAB\cnsp_utils
% addpath ..\..\MATLAB\cnsp_utils\cnd
% addpath ..\..\MATLAB\NoiseTools
% addpath ..\..\MATLAB\eeglab

addpath ./CNSP-resources/CNSP/libs/cnsp_utils
addpath ./CNSP-resources/CNSP/libs/cnsp_utils/cnd
addpath ./CNSP-resources/CNSP/libs/NoiseTools
addpath ./CNSP-resources/CNSP/libs/eeglab

%% Parameters preprocessing
%dataEegFolder = '.\outputs\64Hz\EEG\'; 
%dataCNDFolder = '.\outputs\64Hz\CND\';
% June 19, 2024 changing to new output
%dataEegFolder = './eegProject/datasets/SLdata/';
%dataCNDFolder = './eegProject/datasets/SLdata/1-30Hz/Mastoids/';\
dataEegFolder = './outputs_new/64Hz/EEG/'; 
dataCNDFolder = './outputs_new/64Hz/CND/';
pre = '1-30Hz';
reref_type = 'Mastoids';
downfreq = 64;

%% General parameters
NTRIALS = 14;
conditions = {'V', 'R'};
%subjects = {'698908', '752086'};
subjects = {'164123', '319885', '225637', '917095', '676934', '645737', ...
    '819069', '707776', '558393', '243067', '777647', '727892', '385164', ...
    '613333', '252198', '698908', '888521', '910810',... 
    '992123', '437658', '966706', '780931'};
dataType = 'EEG';
deviceName = 'Biosemi';
chanlocs = load([dataCNDFolder, 'chanlocs64.mat']).chanlocs;

%% Generate data stim for each condition
for condition = conditions
    cond = condition{1};

    % Load data stim for this condition
    load([dataCNDFolder, cond, '/dataStim.mat'], 'stim');
    
    for idxSbj = 1:length(subjects)
        subject = subjects{idxSbj};
        disp(['Processing subject ', subject])
        curr_path = [dataEegFolder,pre,'/',reref_type,'/',subject,'/'];

        % Get preprocessing pipeline from csv
        T = readtable([curr_path, 'preprocessing_pipeline.csv']);
        preprocessingPipeline = {['LPF ', T.LPF, 'Hz'], ... 
                                 ['HPF ', T.HPF, 'Hz'], ...
                                 ['Reref ', T.reref_type], ...
                                 ['Down ', downfreq, 'Hz']};

        % Create data struct for this subject
        eeg = struct();
        eeg.dataType = dataType;
        eeg.deviceName = deviceName;
        eeg.condition = cond;
        eeg.chanlocs = chanlocs;
        eeg.extChan{1,1} = struct('description', 'Mastoids', 'chanlocs', []);
        eeg.fs = downfreq;
        eeg.reRef = T.reref_type;
        eeg.preprocessingPipeline = preprocessingPipeline;
        
        % Fill data structure with trials
        for idxTrial = 1: NTRIALS
            stim_name = stim.audioFiles{idxTrial};
            stim_path = [curr_path,stim_name];
            % Handle missing trials
            if ~exist(stim_path, 'file')
                eeg.nonemptyTrials{idxTrial} = 0; 
                eeg.data{idxTrial} = [];
                eeg.extChan{1,1}.data{idxTrial} = [];
                disp([stim_path, ' does not exist'])
            else
                eeg.nonemptyTrials{idxTrial} = 1;
                load(stim_path);
                eeg.data{idxTrial} = trial_data';
                eeg.extChan{1,1}.data{idxTrial} = trial_mastoids';
            end
        end

        % Interpolate bad channels
        disp('Interpolating bad channels...')
        if isfield(eeg,'chanlocs')
            for i = 1:numel(eeg.data)
                if eeg.nonemptyTrials{i} == 1
                    eeg.data{i} = removeBadChannels(eeg.data{i},eeg.chanlocs);
                end
            end
            eeg = cndNewOp(eeg,'removeBadChannels');
        end

        % Saving preprocessed data
        eegFolder = [dataCNDFolder, cond, '/', pre, '/', reref_type];
        mkdir(eegFolder)
        eegPreFilename = [eegFolder, '/pre_dataSub', num2str(idxSbj), '.mat'];
        disp(['Saving preprocessed EEG data: ', cond, ' of subject ' num2str(idxSbj)])
        save(eegPreFilename,'eeg');

    end

end
disp('Done!')