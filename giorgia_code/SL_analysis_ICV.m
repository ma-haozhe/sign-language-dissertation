%% Init
close all;
clear; clc;

% Set main paths
% base_folder = '.\outputs\30Hz\features\';
% path_fig = '.\outputs\30Hz\figures\'; 
% datafolder = '.\outputs\30Hz\CND\';
% TRFfolder = '.\outputs\30Hz\';

base_folder = './outputs_new/64Hz/features/';
path_fig = './outputs_new/64Hz/figures/';
datafolder = './outputs_new/64Hz/CND/';
TRFfolder = './outputs_new/64Hz/';


%% Parameters
reRefType = 'Mastoids';
cv = 'cvonly';
%June 19, 2024 changed here the lag is -50_450
%lags = '-150_550/';
lags = '-50_450/';
%pre = '0.5-8Hz\';
% what does this pre do here?
pre = '1-30Hz/';

% for viz of the results
conditions = {'V', 'R'}; 
conds_for_plot = {'right', 'reversed'}; 
decoding = 0;
absolute_value_corr = 0;
p_thresh = 0.05;
p_thresh_trfs = 0.05;
do_mafdr = 1;
average_channels = 1;
shuffling = 1;
title_str = 'ICV';
feature_name = 'A2';
feature_name2 = '';
feat_of_interest = 1;  
chan_of_interest = 30;   
chan_of_interest_str = 'Oz';
channels_of_interest = [37, 47, 30];
channels_of_interest_str = {'Fz', 'Cz', 'Pz'};
tmin = -100; 
tmax = 500;
t1 = 50;
t2 = 190;
t3 = 300;
ylims_auto = 1;
ylims = [-1, 2];
common_colorbar = 1;
corr_factor = 1;
set_zlim = 0;
pre_set_zlim_all = 0.06;
pre_set_zlim_dif = 0.01;
avg_over_subject = 0;
avg_over_chs = 1;
brain_area = '';  %centrofrontal-
subjects_to_select = [1,22];  

%% Define parameters figures
xlims = [tmin, tmax];
x = 1:length(conditions);
FontSize = 20;  % 40
FontSizeTitle = 30;  % 55
set(gca,'DefaultTextFontSize',FontSize)
set(0,'defaultfigurecolor',[1 1 1])
alpha = 0.2;
smth = 1;
sem = 1;

% Define RGB values for blue, green, and red and combine them
blue = [0 0.4470 0.7410];
green = [0.4660 0.6740 0.1880];
red = [0.8500 0.3250 0.0980];
colormap_trfs = [blue; green];

% Define the colormap
colormap_topo = flipud(hot(1024));

%% Extract results
rcat = [];
gcat = [];
scat = [];
barv = [];
rcat_dec = [];
gcat_dec = [];
scat_dec = [];
barv_dec = [];
rcat_shu = [];
gcat_shu = [];
scat_shu = [];
barv_shu = [];
rcat_diff = [];
gcat_diff = [];
scat_diff = [];
barv_diff = [];
for idx_cond = 1: length(conditions)
    cond = conditions{idx_cond};
    %% Load EEG data for metadata and channel locs
    %eegPreFilename = [datafolder, cond,'\',pre,reRefType,'\pre_dataSub1.mat'];
    eegPreFilename = [datafolder, cond,'/',pre,reRefType,'/pre_dataSub22.mat'];
    load(eegPreFilename,'eeg')

    %% Load TRFs and correlations
    % Build paths
    tfolder = [TRFfolder, 'TRFs/', pre, lags, reRefType, '/'];
    prename = [feature_name, '_', cond, '_'];
    
    % Load TRFs
    modelAll_path = [tfolder, prename, 'modelAll_', cv, '.mat'];
    load(modelAll_path, 'modelAll');
    model_all{idx_cond} = modelAll;

    % Load Pearson's
    r_all_path = [tfolder, prename, 'rpredAll_', cv, '.mat'];
    load(r_all_path, 'rpredAll');
    r_all{idx_cond} = rpredAll;
    if absolute_value_corr
        r_all{idx_cond} = abs(r_all{idx_cond});
    end

    if decoding
        % Load Pearson's decoding
        tfolder_dec = [TRFfolder, 'Decoders/', pre, lags, reRefType, '/'];
        r_all_path = [tfolder_dec, prename, 'rpredAll_', cv, '.mat'];
        load(r_all_path, 'rpredAll');
        r_dec{idx_cond} = rpredAll;
        if absolute_value_corr
            r_dec{idx_cond} = abs(r_dec{idx_cond});
        end  
    end

    % Load data shuffeling
    if shuffling   
        prename = [feature_name,'_',cond, '_'];
        r_shu_path = [tfolder, prename, 'rpredShu_' cv, '.mat'];
        r_shu{idx_cond} = load(r_shu_path);
        r_shu{idx_cond} = r_shu{idx_cond}.rpredShu;
    else
        prename = [feature_name2, '_', cond, '_'];
        r_shu_path = [tfolder, prename, 'rpredAll_', cv, '.mat'];
        r_shu{idx_cond} = load(r_shu_path, 'rpredAll');
        r_shu{idx_cond} = r_shu{idx_cond}.rpredAll;
    end

    if absolute_value_corr
        r_shu{idx_cond} = abs(r_shu{idx_cond});
    end

    %% Get TRFs for a feature of interest and compute significance
    % Stack normalized models - considering only one feature
    wAvgFeat = [];
    for sub = 1:length(modelAll)
        if isempty(modelAll(sub).w)
            continue
        end
        w = modelAll(sub).w;

        % Select the feature of interest
        w = w(feat_of_interest, :, :);
        
        % normalize
        m = mean(w, 2);
        sd = std(w(:));
        w = w - m;
        w = w/sd; 
        wAvgFeat = [wAvgFeat; w];
    end
    model_avg_feat{idx_cond} = wAvgFeat;

    % ttest on TRF of the feature of interest - variance over subjects
    nlags = size(model_avg_feat{idx_cond}, 2);
    nsbjs = size(model_avg_feat{idx_cond}, 1);
    ttest_vector = [];
    ptest_vector = [];
    if nsbjs > 1
        for lag = 1:nlags
            distr = squeeze(squeeze(model_avg_feat{idx_cond}(:, lag, chan_of_interest)));
            [h,p,ci,stats] = ttest(distr, 0, 'Alpha', 0.05);
            ptest_vector = [ptest_vector; p];
            ttest_vector = [ttest_vector; h];
        end
    end
    ptest_final{idx_cond} = ptest_vector;
    ttest_final{idx_cond} = ttest_vector;

    %% Averaging TRFs and Pearson's across subjects
    % Compute average TRF over subjects
    model_avg_sbj{idx_cond} = mTRFmodelAvg(modelAll,1);

    % Compute average correlation scores across subjects
    rpred_avg_sbj{idx_cond} = squeeze(mean(r_all{idx_cond},1));

    %% Compute statistics on Pearson's for each channel for topographies
    % Get statistics for each channel so to see significant channels
    [tpred{idx_cond}, ppred{idx_cond}] = ttest(r_all{idx_cond}); 

    % Compute average correlation scores (for topographies)
    rdiff{idx_cond} = r_all{idx_cond} - r_shu{idx_cond};
    rdiff_avg_sbj{idx_cond} = squeeze(mean(rdiff{idx_cond},1));

    %% Compute statistics on Pearson's for barplots 
    % (on 1 channel or on the average of all channels)
    if average_channels
        rcat{idx_cond} = squeeze(mean(r_all{idx_cond},2));
    else
        rcat{idx_cond} = r_all{idx_cond}(:, chan_of_interest);
    end 
    [t(idx_cond), p(idx_cond)] = ttest(rcat{idx_cond});
    gcat = [gcat; repmat({conds_for_plot{idx_cond}},length(rcat{idx_cond}),1)];
    scat = [scat; std(rcat{idx_cond})/sqrt(length(rcat{idx_cond}))];
    barv = [barv, mean(rcat{idx_cond})];

    % assign color grey/black according to significnce
    if t(idx_cond) == 0
        C(idx_cond, :) = [.7 .7 .7];
    else
        C(idx_cond, :) = [.0 .0 .0];
    end

    %% Compute statistics on shuffeled Pearson's for barplots 
    % (on 1 channel or on the average of all channels)
    if average_channels
        rcat_shu{idx_cond} = squeeze(mean(r_shu{idx_cond},2));
    else
        rcat_shu{idx_cond} = r_shu{idx_cond}(:, chan_of_interest);
    end 
    [t_shu(idx_cond), p_shu(idx_cond)] = ttest(rcat_shu{idx_cond});
    gcat_shu = [gcat_shu; repmat({conds_for_plot{idx_cond}},length(rcat_shu{idx_cond}),1)];
    scat_shu = [scat_shu; std(rcat_shu{idx_cond})/sqrt(length(rcat_shu{idx_cond}))];
    barv_shu = [barv_shu, mean(rcat_shu{idx_cond})];

    % assign color grey/black according to significnce
    if t_shu(idx_cond) == 0
        C_shu(idx_cond, :) = [.7 .7 .7];
    else
        C_shu(idx_cond, :) = [.0 .0 .0];
    end

    %% Compute statistics on Person's differences with shuffeling
    % (on 1 channel or on the average of all channels)
    if average_channels
        rcat_diff{idx_cond} = mean(r_all{idx_cond},2) - mean(r_shu{idx_cond},2);
    else
        rcat_diff{idx_cond} = r_all{idx_cond}(:, chan_of_interest) - r_shu{idx_cond}(:, chan_of_interest);
    end

    % Get statistics for barplots
    [t_diff(idx_cond), p_diff(idx_cond)] = ttest(rcat_diff{idx_cond});
    gcat_diff = [gcat_diff; repmat({cond},length(rcat_diff{idx_cond}),1)];
    scat_diff = [scat_diff; std(rcat_diff{idx_cond})/sqrt(length(rcat_diff{idx_cond}))];
    barv_diff = [barv_diff, mean(rcat_diff{idx_cond})];  

    % assign color grey/black according to significnce
    if t_diff(idx_cond) == 0
        C_diff(idx_cond, :) = [.7 .7 .7];
    else
        C_diff(idx_cond, :) = [.0 .0 .0];
    end

end

%% Plots topographies and TRFs
fig = figure('Renderer', 'painters', 'Position', [9 9 800 800]); 
n_conditions = length(conditions);
n_rows = 3;

for idx_cond = 1:n_conditions
    %% Person's for each channel on topographies
    subplot(n_rows,n_conditions,idx_cond) 
    title(conds_for_plot{idx_cond}, 'FontSize', FontSize)  

    % Compute significance
    rplot = rpred_avg_sbj{idx_cond};
    [tpred{idx_cond}, ppred{idx_cond}] = ttest(r_all{idx_cond}); 
    if do_mafdr
        pvalues = ppred{idx_cond};
        [~,qvalues] = mafdr(pvalues);
        pvalues_curr = qvalues;
    else
        pvalues_curr = ppred{idx_cond};
    end
    rplot(pvalues_curr>p_thresh) = 0;

    % Get zlim colorbar
    if common_colorbar
        zlim_max = max(cell2mat(rpred_avg_sbj));
    else
        zlim_max = max(rplot);
    end

    if set_zlim
        zlim_max = pre_set_zlim_all;
    end

    % Plot topographies
    topoplot(rplot, eeg.chanlocs,'maplimits',[0, zlim_max],'whitebk','on',...
                'colormap', colormap_topo, 'electrodes', 'off')
    if common_colorbar
        if idx_cond == 1
            c = colorbar;
            c.Location = 'westoutside';
            c.Label.String = ['r'];   
            c.FontSize = FontSize;
        end
    else
        c = colorbar;
        c.Location = 'westoutside';
        c.Label.String = ['r '];   
        c.FontSize = FontSize;
    end

    %% Person's difference for each channel on topographies
    subplot(n_rows,n_conditions,idx_cond+n_conditions)

    % Compute significance
    rplot = rdiff_avg_sbj{idx_cond};
    [tpred{idx_cond}, ppred{idx_cond}] = ttest(rdiff{idx_cond}); 
    if do_mafdr
        pvalues = ppred{idx_cond};
        [~,qvalues] = mafdr(pvalues);
        pvalues_curr = qvalues;
    else
        pvalues_curr = ppred{idx_cond};
    end
    rplot(pvalues_curr>p_thresh) = 0;

    % Get zlim colorbar
    if common_colorbar
        zlim_max = max(cell2mat(rdiff_avg_sbj));
    else
        zlim_max = max(rplot);
    end

    if set_zlim
        zlim_max = pre_set_zlim_dif;
    end

    % Plot topographies
    topoplot(rplot, eeg.chanlocs,'maplimits',[0, zlim_max],'whitebk','on',...
                'colormap', colormap_topo, 'electrodes', 'off')
    if common_colorbar
        if idx_cond == 1
            c = colorbar;
            c.Location = 'westoutside';
            c.Label.String = ['r gain'];   
            c.FontSize = FontSize;
        end
    else
        c = colorbar;
        c.Location = 'westoutside';
        c.Label.String = ['r gain'];   
        c.FontSize = FontSize;
    end

    %% Plot TRFs for a given feature and channel with significance
    subplot(n_rows,n_conditions,idx_cond+(2*n_conditions))
    ax = stdshade(squeeze(model_avg_feat{idx_cond}(:, :, chan_of_interest)), ...
                  alpha,colormap_trfs(idx_cond, :),model_avg_sbj{idx_cond}.t,smth,sem);
    if ylims_auto
        ylim auto
    else
        ylim(ylims) 
    end
    hold on
    ttest_final{idx_cond}(ptest_final{idx_cond} > p_thresh_trfs) = 0;
    ttest_trfs = ttest_final{idx_cond};
    shaded_patch_significant_timepoints(model_avg_sbj{idx_cond}.t, ttest_trfs, ax)
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5) 
    xlabel('Time lag (ms)', 'FontSize', FontSize) 
    if idx_cond == 1
    	ylabel(['TRF at ', chan_of_interest_str], 'FontSize', FontSize)
    end
    xlim(xlims)
    xticks([0 200 400 600])
    xticklabels({'0','0.2','0.4', '0.6'})
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',FontSize)
    set(gca,'XTickLabelMode','auto')
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',FontSize)
    set(gca,'YTickLabelMode','auto')
end

set(gcf,'color','w');


%% Plots topographies and TRFs for each subject
for sbj_id = 1:nsbjs
    fig = figure('Renderer', 'painters', 'Position', [9 9 700 700]); 
    n_conditions = length(conditions);
    n_rows = 3;

    for idx_cond = 1:n_conditions
        %% Person's for each channel on topographies
        subplot(n_rows,n_conditions,idx_cond) 
        title(conds_for_plot{idx_cond}, 'FontSize', FontSize)  
    
        rplot = r_all{idx_cond}(sbj_id, :);
        zlim_max = max(rplot);
        topoplot(rplot, eeg.chanlocs,'maplimits',[0, zlim_max],'whitebk','on','colormap', colormap_topo, 'electrodes', 'off')
        c = colorbar;
        c.Location = 'westoutside';
        c.Label.String = ['r '];   
        c.FontSize = FontSize;
    
        %% Person's difference for each channel on topographies
        subplot(n_rows,n_conditions,idx_cond+n_conditions)

        rplot = r_all{idx_cond}(sbj_id, :)-r_shu{idx_cond}(sbj_id, :);
        zlim_max = max(rplot);
        topoplot(rplot, eeg.chanlocs,'maplimits',[0, zlim_max],'whitebk','on','colormap',colormap_topo, 'electrodes', 'off')
        c = colorbar;
        c.Location = 'westoutside';
        c.Label.String = ['r gain'];   
        c.FontSize = FontSize;

        %% Plot TRFs for a given feature and channel with significance
        subplot(n_rows,n_conditions,idx_cond+(2*n_conditions))
        curr_model = model_all{idx_cond}(sbj_id).w;
        m = mean(curr_model, 2);
        curr_model = curr_model - m;
        sd = std(curr_model(:));
        curr_model = curr_model/sd; 
        ax = stdshade(squeeze(curr_model(feat_of_interest, :, chan_of_interest)), ...
                      alpha,colormap_trfs(idx_cond, :),model_avg_sbj{idx_cond}.t,smth,sem);
        if ylims_auto
            ylim auto
        else
            ylim(ylims) 
        end
        hold on
        ttest_final{idx_cond}(ptest_final{idx_cond} > p_thresh_trfs) = 0;
        ttest_trfs = ttest_final{idx_cond};
        shaded_patch_significant_timepoints(model_avg_sbj{idx_cond}.t, ttest_trfs, ax)
        yline(0, '-', 'Alpha', 0.5)
        xline(0, '-', 'Alpha', 0.5) 
        xlabel('Time lag (ms)', 'FontSize', FontSize) 
        if idx_cond == 1
    	    ylabel(['TRF at ', chan_of_interest_str], 'FontSize', FontSize)
        end
        xlim(xlims)
        xticks([0 200 400 600])
        xticklabels({'0','0.2','0.4', '0.6'})
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',FontSize)
        set(gca,'XTickLabelMode','auto')
        a = get(gca,'YTickLabel');  
        set(gca,'YTickLabel',a,'fontsize',FontSize)
        set(gca,'YTickLabelMode','auto')
    end
end

set(gcf,'color','w');

%% Plot barplots with statistics
fig = figure('Renderer', 'painters', 'Position', [10 10 800 800]); 
% sgtitle(title_str, 'FontSize', FontSizeTitle)

bar = [barv_shu(1); barv(1); 0; barv_shu(2); barv(2)];
sca = [scat_shu(1); scat(1); 0; scat_shu(2); scat(2)];
col = [[C_shu(1, :)]; [C(1, :)]; [1, 1, 1]; [C_shu(2, :)]; [C(2, :)]];
xlabel_bars = {'shuff', 'IVC', ' ', 'shuff', 'IVC'};

% Compute cross-significance for Pearson's
P = nan(5, 5);
[h_tmp, P(1,2)] = ttest(cell2mat(rcat(1)), cell2mat(rcat_shu(1)));
[h_tmp, P(4,5)] = ttest(cell2mat(rcat(2)), cell2mat(rcat_shu(2)));
% [h_tmp, P(1,4)] = ttest(cell2mat(rcat_shu(2)), cell2mat(rcat_shu(1)));
% [h_tmp, P(2,5)] = ttest(cell2mat(rcat(2)), cell2mat(rcat(1)));
P = P.*corr_factor;
% P(1,2) =  p_diff(idx_cond)
PT = P'; 
lidx = tril(true(size(P)), -1);
P(lidx) = PT(lidx);

% Plot barplot Pearson's
superbar(bar, 'E', sca, 'P', P, 'BarFaceColor', col, ...
    'Orientation', 'v', 'ErrorbarStyle', 'T', 'PLineOffset', 0.005);
hold on
x_coo_shu = ones(size(rcat_shu{1})).*(1+(rand(size(rcat_shu{1}))-0.5)/10);
x_coo_dis = ones(size(rcat{1})).*(2+(rand(size(rcat{1}))-0.5)/10);
y_coo_shu = rcat_shu{1};
y_coo_dis = rcat{1};
for s = 1:length(subjects_to_select)
    x = [x_coo_shu(s); x_coo_dis(s)];
    y = [y_coo_shu(s); y_coo_dis(s)];
    plot(x, y,'r')
%     plot(x, y, 'DisplayName', num2str(s))
end
scatter(x_coo_shu,y_coo_shu,'r','filled')
scatter(x_coo_dis,y_coo_dis,'r','filled')

x_coo_shu = ones(size(rcat_shu{2})).*(4+(rand(size(rcat_shu{2}))-0.5)/10);
x_coo_dis = ones(size(rcat{2})).*(5+(rand(size(rcat{2}))-0.5)/10);
y_coo_shu = rcat_shu{2};
y_coo_dis = rcat{2};
for s = 1:length(subjects_to_select)
    x = [x_coo_shu(s); x_coo_dis(s)];
    y = [y_coo_shu(s); y_coo_dis(s)];
    plot(x, y,'r')
%     plot(x, y, 'DisplayName', num2str(s))
end
scatter(x_coo_shu,y_coo_shu,'r','filled')
scatter(x_coo_dis,y_coo_dis,'r','filled')

xlabel('      V               R', 'FontSize', FontSize)
ylabel('r', 'FontSize', FontSize) 
ylim auto
a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',FontSize)
set(gca,'YTickLabelMode','auto')
xtick = 1:5;
set(gca,'xtick',xtick,'xticklabel',xlabel_bars, 'FontSize', FontSize)
set(gcf,'color','w');
% legend({'p<0.05'})
% legend('boxoff')

%% Plot TRFs and weights on topographies
fig = figure();  %'Renderer', 'painters', 'Position', [10 10 350 800]);
sgtitle(title_str)

n_rows = 4;
for idx_cond = 1:n_conditions
    subplot(n_rows,n_conditions,idx_cond) 

    % Plot TRF
    avgModel = model_avg_sbj{idx_cond};
    mTRFplot(avgModel, 'trf', feat_of_interest, chan_of_interest, xlims);
    title(conds_for_plot{idx_cond}, 'FontSize', FontSize)
    xlabel('Time lag (ms)', 'FontSize', FontSize)
    ylabel('amplitude (a.u.)', 'FontSize', FontSize)
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5) 
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    xlim(xlims)
    ylim auto
%     axis square
    
    % xlim
    [topo_minValue,idx_min] = min(abs(avgModel.t-tmin));   
    [topo_maxValue,idx_max] = min(abs(avgModel.t-tmax));

    % ylim
    ylim_max = max(max(abs(avgModel.w(:,idx_min:idx_max,:)),[],3),[],2);
    ylim_min = min(min(avgModel.w(:,idx_min:idx_max,:),[],3),[],2);
    
    % zlim
    lim = max(max(abs(avgModel.w(feat_of_interest,idx_min:idx_max,:)),[],3),[],2);
    
    % topo time values
    [topo1_val, topo1_idx] = min(abs(avgModel.t - t1));
    [topo2_val, topo2_idx] = min(abs(avgModel.t - t2));
    [topo3_val, topo3_idx] = min(abs(avgModel.t - t3));
    
    % Plot avg TRF model
    subplot(n_rows,n_conditions,idx_cond+n_conditions)
    topoplot(avgModel.w(feat_of_interest,topo1_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(avgModel.t(topo1_idx)),' ms'], 'FontSize', FontSize)

    subplot(n_rows,n_conditions,idx_cond+(2*n_conditions))
    topoplot(avgModel.w(feat_of_interest,topo2_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(avgModel.t(topo2_idx)),' ms'], 'FontSize', FontSize)
    
    subplot(n_rows,n_conditions,idx_cond+(3*n_conditions))
    topoplot(avgModel.w(feat_of_interest,topo3_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(avgModel.t(topo3_idx)),' ms'], 'FontSize', FontSize)
end
% print(fig, [path_fig, feature_name, '_TRFs_'],'-r800','-dpng');


%% Plot TRFs and weights on topographie - one plot for each condition
n_rows = 3;
for idx_cond = 1:n_conditions
    figure('Renderer', 'painters')
    subplot(n_rows,3,[1:6]) 
    avgModel = model_avg_sbj{idx_cond};
    ax = plot(avgModel.t,squeeze(avgModel.w(feat_of_interest,:,:)));
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5)  
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    title(conds_for_plot{idx_cond}, 'FontSize', FontSizeTitle)
    xlabel('time lag (ms)', 'FontSize', FontSize)
    if idx_cond == 1
        ylabel('TRF amplitude (a.u.)', 'FontSize', FontSize)
    end
    xlim(xlims)
    if ylims_auto
        ylim auto
    else
        ylim(ylims) 
    end
    hold on
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',FontSize)
    set(gca,'XTickLabelMode','auto')
    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',FontSize)
    set(gca,'YTickLabelMode','auto')
    set(gcf,'color','w');
 
    % xlim
    [topo_minValue,idx_min] = min(abs(avgModel.t-tmin));   
    [topo_maxValue,idx_max] = min(abs(avgModel.t-tmax));

    % ylim
    ylim_max = max(max(abs(avgModel.w(:,idx_min:idx_max,:)),[],3),[],2);
    ylim_min = min(min(avgModel.w(:,idx_min:idx_max,:),[],3),[],2);
    
    % zlim
    lim = max(max(abs(avgModel.w(feat_of_interest,idx_min:idx_max,:)),[],3),[],2);
    
    % topo time values
    [topo1_val, topo1_idx] = min(abs(avgModel.t - t1));
    [topo2_val, topo2_idx] = min(abs(avgModel.t - t2));
    [topo3_val, topo3_idx] = min(abs(avgModel.t - t3));
    
    % Plot avg TRF model
    subplot(n_rows,3,7)
    topoplot(avgModel.w(feat_of_interest,topo1_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(t1),' ms'], 'FontSize', FontSize, 'Units', 'normalized', 'Position', [0.5, -0.25, 0]);
    c = colorbar;
    c.Location = 'westoutside';
    c.Label.String = ['weights'];   
    c.FontSize = FontSize;

    subplot(n_rows,3,8)
    topoplot(avgModel.w(feat_of_interest,topo2_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(t2),' ms'], 'FontSize', FontSize, 'Units', 'normalized', 'Position', [0.5, -0.25, 0]);
    
    subplot(n_rows,3,9)
    topoplot(avgModel.w(feat_of_interest,topo3_idx,:),eeg.chanlocs,'maplimits',[-lim,lim],'whitebk','on')
    title([num2str(t3),' ms'], 'FontSize', FontSize, 'Units', 'normalized', 'Position', [0.5, -0.25, 0]);

    set(gcf,'color','w');
end

%% Plot butterfly plots
fig = figure();  %'Renderer', 'painters', 'Position', [10 10 350 800]);
sgtitle(title_str)

n_rows = 2;
for idx_cond = 1:n_conditions
    subplot(n_rows,n_conditions,idx_cond) 
    avgModel = model_avg_sbj{idx_cond};
    plot(avgModel.t,squeeze(avgModel.w(feat_of_interest,:,:)))
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5)  
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    title(['TRF ', conds_for_plot{idx_cond}], 'FontSize', FontSize)
    xlabel('Time lag (ms)', 'FontSize', FontSize)
    ylabel('amplitude (a.u.)', 'FontSize', FontSize)
    xlim(xlims)
    ylim auto

    subplot(n_rows,n_conditions,idx_cond+n_conditions) 
    plot(avgModel.t,std(avgModel.w(feat_of_interest,:,:),[],3))
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5)  
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    title(['GFP ', conds_for_plot{idx_cond}], 'FontSize', FontSize)
    xlabel('Time lag (ms)', 'FontSize', FontSize)
    ylabel('amplitude (a.u.)', 'FontSize', FontSize)
    xlim(xlims)
    ylim auto
end

%% Plot TRF of channels of interest
fig = figure();  %'Renderer', 'painters', 'Position', [10 10 350 800]);
tiledlayout(1, length(channels_of_interest), 'TileSpacing','tight');
% sgtitle(title_str)

n_rows = 1;
for idx_ch = 1:length(channels_of_interest)
    nexttile
    for idx_cond = 1:n_conditions    
        stdshade(squeeze(model_avg_feat{idx_cond}(:, :, channels_of_interest(idx_ch))), ...
                 alpha,colormap_trfs(idx_cond, :),model_avg_sbj{idx_cond}.t,smth,sem);
        hold on
    end
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5) 
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    legend('speech', '', 'sung speech', '', '')
    title(channels_of_interest_str(idx_ch), 'FontSize', FontSize)
    xlabel('Time lag (ms)', 'FontSize', FontSize)
    ylabel('amplitude (a.u.)', 'FontSize', FontSize)
    xlim(xlims)
    ylim auto
    axis square
end

%% Plot TRF of channels of interest with Cohen's d
n_rows = 3;
for idx_ch = 1:length(channels_of_interest)
    fig = figure('Renderer', 'painters');
    subplot(n_rows,1,[1:2]) 
    for idx_cond = 1:n_conditions    
        stdshade(squeeze(model_avg_feat{idx_cond}(:, :, channels_of_interest(idx_ch))), ...
                 alpha,colormap_trfs(idx_cond, :),model_avg_sbj{idx_cond}.t,smth,sem);
        hold on
    end
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5) 
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 
    if idx_ch == 3
        legend({'speech', '', 'sung speech', '', ''}, 'FontSize', FontSize)
    end
    title(channels_of_interest_str(idx_ch), 'FontSize', FontSize)
    if idx_ch == 1
        ylabel('TRF amplitude (a.u.)', 'FontSize', FontSize)
        yticks([-2 -1 0 1 2]);
        a = get(gca,'YTickLabel');  
        set(gca,'YTickLabel',a,'fontsize',FontSize)
    else
        yticklabels({});
    end
    xticklabels({});
    xlim(xlims)
    ylim(ylims) 

    % Compute effect size
    d_song_speech = [];
    if size(model_avg_feat{1}, 3) > 1
        for lag = 1:size(model_avg_feat{1}, 2)
            distr_song   = squeeze(model_avg_feat{1}(:, lag, channels_of_interest(idx_ch)));
            distr_speech = squeeze(model_avg_feat{2}(:, lag, channels_of_interest(idx_ch)));
            [h, p, ci, stats] = ttest2(distr_song, distr_speech, 'vartype', 'unequal');
            d_song_speech = [d_song_speech; abs(stats.tstat / sqrt(stats.df + 1))];
        end
    end

    subplot(n_rows,1,3) 
    plot(avgModel.t,d_song_speech, 'color', blue, 'LineWidth', 2);
    hold on
    xlabel('Time lag (ms)', 'FontSize', FontSize)
    if idx_ch == 1
        ylabel('Effect size', 'FontSize', FontSize)
        yticks([0.2 0.5 0.8]);
        a = get(gca,'YTickLabel');  
        set(gca,'YTickLabel',a,'fontsize',FontSize)
    else
        yticklabels({});
    end
    xlim(xlims)
    ylim([0, 1])
    hold on
    yline(0, '-', 'Alpha', 0.5)
    xline(0, '-', 'Alpha', 0.5) 
    xline(t1, '--', 'Alpha', 0.9) 
    xline(t2, '--', 'Alpha', 0.9) 
    xline(t3, '--', 'Alpha', 0.9) 

    yline(0.2, '-', 'LineWidth', 2)  % 'Small effect size'
    yline(0.5, '-', 'LineWidth', 2)  % 'Medium effect size'
    yline(0.8, '-', 'LineWidth', 2)  % 'Large effect size'

    if idx_ch == 3
        legend({'song-speech'}, 'FontSize', FontSize)
    end

    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',FontSize)
    set(gca,'XTickLabelMode','auto')
    set(gcf,'color','w');
end

%% Plots difference topographies - correlations
fig = figure(); 

rplot_song   = rpred_avg_sbj{2};
rplot_speech = rpred_avg_sbj{1};
rdiff_song_speech = rplot_song - rplot_speech;

rall_song   = r_all{2};
rall_speech = r_all{1};
[tpred_song_speech, ppred_song_speech] = ttest(rall_song-rall_speech); 

% ppred_song_melody = ppred_song_melody*64;
% ppred_song_speech = ppred_song_speech*64;
% [~,ppred_song_melody] = mafdr(ppred_song_melody);
% [~,ppred_song_speech] = mafdr(ppred_song_speech);

rdiff_song_speech(ppred_song_speech>p_thresh) = 0;

r_diff{2} = rdiff_song_speech;

% Get zlim colorbar
zlim_max = max(max(cell2mat(r_diff)));
zlim_min = min(min(cell2mat(r_diff)));

subplot(1,2,2) 
title('Song-speech', 'FontSize', FontSizeTitle)
rplot = rdiff_song_speech;
if common_colorbar
    if zlim_max<=zlim_min
        zlim_min = 0;
        zlim_max = 1;
    end
    lims = [zlim_min, zlim_max];
else
    zlim_min = min(rplot);
    zlim_max = max(rplot);
    if zlim_max<=zlim_min
        zlim_min = 0;
        zlim_max = 1;
    end
    lims = [zlim_min, zlim_max];
end
topoplot(rplot, eeg.chanlocs,'maplimits',lims,'whitebk','on')
c = colorbar;
set(c, 'Orientation', 'horizontal');
c.Location = 'westoutside';
c.Label.String = ['Pearson correlation (r)'];   
c.FontSize = FontSize;

set(gcf,'color','w');


%% Plots difference topographies - TRF weights
figure('Renderer', 'painters')

[topo1_val, topo1_idx] = min(abs(avgModel.t - t1));
[topo2_val, topo2_idx] = min(abs(avgModel.t - t2));
[topo3_val, topo3_idx] = min(abs(avgModel.t - t3));
topo_idxs = [topo1_idx, topo2_idx, topo3_idx];
time_values = [t1, t2, t3];

TRF_song   = model_avg_sbj{1}.w;
TRF_speech = model_avg_sbj{2}.w;
TRFdiff_song_speech = TRF_song - TRF_speech;

i = 1;
for lag = topo_idxs
    p_song_melody = [];
    p_song_speech = [];
    for ch = 1:size(model_avg_feat{1}, 3)
        distr_song   = squeeze(model_avg_feat{1}(:, lag, ch));
        distr_speech = squeeze(model_avg_feat{2}(:, lag, ch));

        [h, p, ci, stats] = ttest2(distr_song, distr_speech, 'Alpha', 0.05);
        p_song_speech = [p_song_speech; p];
    end
    [~, p_s_s{i}] = mafdr(p_song_speech);
%     p_s_m{i} = p_song_melody*64;
%     p_s_s{i} = p_song_speech*64;
    i = i + 1;
end

for idx_lag = 1:3
    subplot(1,3,idx_lag)
    % Calculate FDR adjusted p-values using the mafdr function
    TRFtoviz = abs(TRFdiff_song_speech(feat_of_interest,topo_idxs(idx_lag),:));
    TRFtoviz(p_s_s{idx_lag}>p_thresh) = 0;
    if common_colorbar
        TRF_diff{1} = TRFdiff_song_speech(feat_of_interest,topo_idxs,:);
        zlim_max = max(max(abs(cell2mat(TRF_diff))));
        zlim_min = min(min(abs(cell2mat(TRF_diff))));
    else
        zlim_max = max(abs(TRFtoviz));
        zlim_min = min(abs(TRFtoviz));
    end
    if zlim_max<=zlim_min
        zlim_min = 0;
        zlim_max = 1;

    end
    topoplot(TRFtoviz,eeg.chanlocs,'maplimits',[zlim_min,zlim_max],'whitebk','on')
    title([num2str(time_values(idx_lag)), ' ms'], 'FontSize', FontSizeTitle) 
    c = colorbar;
    set(c, 'Orientation', 'horizontal');
    c.Location = 'westoutside';
    c.Label.String = ['weights'];   
    c.FontSize = FontSize;
end

set(gcf,'color','w');

%% Topografia
figure;
topoplot(ones(64,1),eeg.chanlocs,'electrodes','numbers')