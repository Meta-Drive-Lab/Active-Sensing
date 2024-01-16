clear;

rootdir = 'C:/Users/A/Documents/MATLAB/driving_model/ActiveSensing';
addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
set(0,'defaultAxesFontSize',20);



% Initialize variables
excludeLastFix = 1;   % exclude last fixation duration, since it's cut short by choosing
useBinnedX = 1;   % Create time bins for switch prob data for smoother results
useBonferroni = 1;  

% Load behavioral data
% taskstr = 'aVBDM';
load(fullfile(data_dir,'datastruct_hkust0902double_change.mat'));

dt = 0.001;
d = 0.0002;
sig2_origModel = 0.02^2;
aGamma_k = 0.3;
k = d/dt;
sig2_k = sig2_origModel/dt;
decbound = 1;


% Get empty datastruct for model
numsubs = length(dstruct_real);
z_reps = 40;  
z_all = repmat(permn(1:3,2),z_reps,1); 
dstruct_addm = getFilledDstruct(dstruct_real,numsubs,z_all);


% Get empirical distribution of fixation behavior for each difficulty level
fixdist = getEmpiricalFixationDist(dstruct_real);
maxdectime = 8;

% RT bin info
binstep = 0.5;
binedges_rt = 1:binstep:4;
binctr_rt = binedges_rt(1:end-1)+binstep/2;

% Get all unique value difference
valdiff_unique = unique(getValdiff(dstruct_real))';

% Simulate behavior
every10Prc = floor(prctile(1:numsubs,10:10:100));
for s = 1:numsubs
    % Generate trials - draw from the prior distribution from the optimal model
    numtrials = length(dstruct_addm(s).trialnum);
    for ti = 1:numtrials
        % aDDM: simulate behavior for this decision boundary
        thisvaldiff = abs(dstruct_addm(s).itemval(ti,1)-dstruct_addm(s).itemval(ti,2));
        if thisvaldiff > 2, thisvaldiff = 2; end
        fixdist_first = fixdist.all_first{fixdist.valdiff==round(thisvaldiff)};
        fixdist_mid = fixdist.all_mid{fixdist.valdiff==round(thisvaldiff)};
        if rand <= 0.5, y = 1;
        else, y=2;
        end
        % sequence of fixations
        yseq_toMaxTime = nan(1,maxdectime/dt);
        ni = 1;
        while sum(~isnan(yseq_toMaxTime)) < maxdectime/dt
            if ni==1, itemfixdur = fixdist_first(randperm(length(fixdist_first),1));
            else
                itemfixdur = fixdist_mid(randperm(length(fixdist_mid),1));
            end
            itemfixN = round(itemfixdur/dt);
            yseq_toMaxTime(ni:ni+itemfixN-1) = y;
            ni = ni+itemfixN;
            y = 3-y;
        end
        yseq_toMaxTime = yseq_toMaxTime(1:maxdectime/dt);
        % Re-write choice
        [dstruct_addm(s).choice(ti),dstruct_addm(s).rt(ti),dstruct_addm(s).fixitem{ti},dstruct_addm(s).fixdur{ti},dstruct_addm(s).tItem(ti,:),~] = run_aDDM(dstruct_addm(s).itemval(ti,:),aGamma_k,dt,k,sig2_k,decbound,yseq_toMaxTime);
    end
end


% Combine all datastructs to plot them efficiently
allDstructs = {dstruct_real,dstruct_addm};
titles_dstruct = {'Data','aDDM'};
colors_dstruct = {'m','b'};
marksize = 10;

% Load params
params.dt = 0.05
RTcutoff = 5;  % Cutoff of RT for trials
switchBinaryMat = struct;  % binary matrix with 1's at time point when switch ocurred
switchBinaryMat.all = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valdiff_low = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.valdiff_hi = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
switchBinaryMat.logreg = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
% Fixation duration - store the fixation duration, stored at fixation onset
% time
fixdurMat = zeros(numsubs,RTcutoff/params.dt,length(allDstructs));
fixdur_all = {};  % Accumulate all middle fixation durations
% Switch number (normalized) for valdiff 
switchnum_valdiff = {}; 
% Switch rate (#switches/rt) for valdiff 
switchrate_valdiff = {}; 
% x-axis
xaxis_orig = 0:params.dt:RTcutoff-params.dt;

% Get switch probability and fixation duration across time
for d = 1:length(allDstructs)
    dstruct = allDstructs{d};
    fixdur_all{d} = [];
    % Switch proportion for value diff & value sum
    itemval_all = [];
    for s = 1:length(dstruct)
        itemval_all = cat(1,itemval_all,dstruct(s).itemval);
    end
    valdiff_unique = unique(abs(itemval_all(:,1)-itemval_all(:,2)))';
    switchnum_valdiff_d = nan(length(dstruct),length(valdiff_unique));
    switchrate_valdiff_d = nan(length(dstruct),length(valdiff_unique));
    for s = 1:length(dstruct)
        % All trials
        numtrials = length(dstruct(s).fixdur);
        switchBinaryMat_sub = zeros(numtrials,RTcutoff/params.dt);
        fixdurMat_sub = nan(numtrials,RTcutoff/params.dt);
        for t = 1:length(dstruct(s).fixdur)
            fixdur_trial = dstruct(s).fixdur{t};
            % Switch probability
            switchTimeIndex = round(round(cumsum(fixdur_trial),2)/params.dt) + 1;  % +1 to include RT of 0
            switchTimeIndex(switchTimeIndex>RTcutoff/params.dt) = [];  % remove any time points greater than maxdectime
            if excludeLastFix == 1 && ~isempty(switchTimeIndex)
                switchTimeIndex(end) = [];
                fixdur_trial(end) = [];
            end
            if ~isempty(switchTimeIndex)
                switchBinaryMat_sub(t,switchTimeIndex) = 1;
            end
            
            % Fixation duration
            fixonsetTimeIndex = [1,switchTimeIndex(1:end-1)];
            if length(fixonsetTimeIndex) > 1
                fixdur_this = fixdur_trial(1:length(switchTimeIndex));
                fixdur_this(1) = [];  % remove first fixation, just look at middle fixations
                fixonsetTimeIndex(1) = [];
                fixdurMat_sub(t,fixonsetTimeIndex) = fixdur_this;
            end
        end
        switchBinaryMat.all(s,:,d) = mean(switchBinaryMat_sub,1);
        % Split trials by upper third and lower third in terms of valdiff
        valdiff = abs(dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2));
        i_lower34 = valdiff < prctile(valdiff,34);
        i_middle50 = valdiff >= prctile(valdiff,34) & valdiff <= prctile(valdiff,66);
        i_upper66 = valdiff > prctile(valdiff,66);
        switchBinaryMat.valdiff_low(s,:,d) = mean(switchBinaryMat_sub(i_lower34,:),1);
        switchBinaryMat.valdiff_mid(s,:,d) = mean(switchBinaryMat_sub(i_middle50,:),1);
        switchBinaryMat.valdiff_hi(s,:,d) = mean(switchBinaryMat_sub(i_upper66,:),1);
        
                
        % Fixation duration
        fixdurMat(s,:,d) = nanmean(fixdurMat_sub,1);
        fixdur_all{d} = cat(1,fixdur_all{d},fixdurMat_sub(~isnan(fixdurMat_sub)));
        
        % 2. Switch proportion vs valdiff
        % Switch rate: #switches / RT
        switchnum = cellfun('length',dstruct(s).fixdur) - 1;
        valdiff_sub = abs(dstruct(s).itemval(:,1)-dstruct(s).itemval(:,2));
        for vd = valdiff_unique
            i_vd = valdiff_sub==vd;
            switchnum_valdiff_d(s,vd == valdiff_unique) = nanmean(switchnum(i_vd)) / mean(switchnum);  % Normalize by mean number of all switches for subject
            rt_without_lastfix = dstruct(s).rt(i_vd) - cellfun(@(v) v(end), dstruct(s).fixdur(i_vd));
            switchrate_valdiff_d(s,vd == valdiff_unique) = nanmean(switchnum(i_vd)./rt_without_lastfix);
        end
        
    end
    switchnum_valdiff{d} = switchnum_valdiff_d;
    switchrate_valdiff{d} = switchrate_valdiff_d;
end

% Create time bins for switch prob data for smoother results
if useBinnedX == 1
    switchBinaryMat_orig = switchBinaryMat;
    binSize = 4;
    i_bin = discretize(xaxis_orig,length(xaxis_orig)/binSize);
    % Get new x-axis, using bin mean
    xaxis_binned = nan(1,length(unique(i_bin)));
    for b = 1:max(i_bin)
        xaxis_binned(b) = mean(xaxis_orig(b==i_bin));
    end
    % Bin data
    for ff = fieldnames(switchBinaryMat)'
        tempMat = nan(size(switchBinaryMat.(char(ff))(:,:,d),1),max(i_bin),size(switchBinaryMat.(char(ff)),3));
        for d = 1:length(allDstructs)
            for b = 1:max(i_bin)
                tempMat(:,b,d) = mean(switchBinaryMat.(char(ff))(:,b==i_bin,d),2);
            end
        end
        switchBinaryMat.(char(ff)) = tempMat;
    end
    xaxis = xaxis_binned;
else
    xaxis = xaxis_orig;
end

% Plot results
% Switch probability
ylimit_switchprob = [0,0.3]; 
xlimit_switchprob = [0,RTcutoff];
for d = 1:length(allDstructs)
    fh = figure('units','normalized','outerposition',[0,0,1,0.4]);
    % All trials
    subplot(1,2,1); pbaspect([1,1,1]); hold on;
    [data_mean,data_se] = getMeanAndSE(switchBinaryMat.all(:,:,d));
    errorbar(xaxis,data_mean,data_se,'Color',colors_dstruct{d});
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title(titles_dstruct{d});
    
    
    % Split by valdiff
    subplot(1,2,2); pbaspect([1,1,1]); hold on;
    [mean_low,se_low] = getMeanAndSE(switchBinaryMat.valdiff_low(:,:,d));
    ax_low = errorbar(xaxis,mean_low,se_low,'g');
    [mean_hi,se_hi] = getMeanAndSE(switchBinaryMat.valdiff_hi(:,:,d));
    ax_hi = errorbar(xaxis,mean_hi,se_hi,'r');
    set(gca,'ylim',ylimit_switchprob,'xlim',xlimit_switchprob);
    xlabel('Time (s)'); ylabel('Switch probability');
    title('Split by ease level');
    % Stats
    [~,pvals] = ttest(switchBinaryMat.valdiff_low(:,:,d),switchBinaryMat.valdiff_hi(:,:,d));
    if useBonferroni==1
        i_sig = pvals < 0.05/size(switchBinaryMat.all,1);
    else
        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);  % FDR-corrected for multiple comparisons
        i_sig = adj_p < 0.05;
    end
    % y_asterisk = 0.14;  % Where to show stars for significance
    % plot(xaxis(i_sig),repmat(y_asterisk,1,sum(i_sig)),'k*');
    legend([ax_low,ax_hi],{'Desicion Ease Level = 0','Desicion Ease Level = 2'});

      
    
end

