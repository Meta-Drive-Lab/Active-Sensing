clear;

rootdir = 'C:/Users/A/Documents/MATLAB/driving_model/ActiveSensing';
addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
set(0,'defaultAxesFontSize',10);


% Set aDDM parameters
dt = 0.001;
d = 0.0002;
sig2_origModel = 0.02^2;
aGamma_k = 0.3;
k = d/dt;
sig2_k = sig2_origModel/dt;
decbound = 1;

% Plotting options
plotOptions = struct;
plotOptions.saveFigures = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behavioral plots for human behavior & aDDM simulations
% Load behavioral data
load(fullfile(data_dir,sprintf('datastruct_hkust0902double_change.mat')));

% Generate behavior using model
numsubs = length(dstruct_real);

whichModels = [0, 1]; % 0: Human data, 1: aDDM simulation
for i = 1:length(whichModels)
   whichModel = whichModels(i);
   % whichModel = 0; % 0: Human data, 1: aDDM simulation
    if whichModel == 0  % Human data
        dstruct = dstruct_real;
        % Plot behavioral plots
        bout_sim = getbehavoutput(dstruct);
        stats_behav = makeBehavPlots_realData(bout_sim,plotOptions);
        actual_data = bout_sim.item1chosen_valdiff_all;

    elseif whichModel == 1  % aDDM model    
        % Get empty datastruct for models
        z_reps = 40;  % Number of iteration
        z_all = repmat(permn(1:3,2),z_reps,1); %Generate all possible binary groups from 1 to 3
        dstruct = getFilledDstruct(dstruct_real,numsubs,z_all);
        count = 0;
        % Simulate behavior
        for s = 1:numsubs
            % Get empirical distribution of fixation behavior for each difficulty level
            fixdist = getEmpiricalFixationDist(dstruct_real);
            maxdectime = 8;
            numtrials = length(dstruct(s).trialnum);
            for ti = 1:numtrials
                thisvaldiff = abs(dstruct(s).itemval(ti,1)-dstruct(s).itemval(ti,2));
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
                
                [dstruct(s).choice(ti),dstruct(s).rt(ti),dstruct(s).fixitem{ti},dstruct(s).fixdur{ti},dstruct(s).tItem(ti,:),~] = run_aDDM(dstruct(s).itemval(ti,:),aGamma_k,dt,k,sig2_k,decbound,yseq_toMaxTime);
                
            end
        end
    bout_sim = getbehavoutput(dstruct);
    stats_behav = makeBehavPlots_model(bout_sim,plotOptions);
    stimulus_data = bout_sim.item1chosen_valdiff_all;
    end

end

% Calculate MSE between simulated and real data
mse_item1chosen = mse(stimulus_data - actual_data);
disp(['MSE between simulated and real data: ', num2str(mse_item1chosen)]);