clear;

rootdir = 'C:/Users/A/Documents/MATLAB/driving_model/ActiveSensing';
addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
set(0,'defaultAxesFontSize',10);


% Set aDDM parameters
dt = 0.001;
d = 0.0002;
sig2_origModel = 0.02^2;
aGamma_k = 0.3;%the adopted parameters are d=0.002, σ=0.02, and θ=0.3
k = d/dt;
sig2_k = sig2_origModel/dt;
decbound = 1;

% Plotting options
plotOptions = struct;
plotOptions.saveFigures = 0;


%% Behavioral plots for human behavior & aDDM simulations
% Load behavioral data
load(fullfile(data_dir,sprintf('datastruct_hkust0902double_change.mat')));

% Generate behavior using model
numsubs = length(dstruct_real);
 
% Get empty datastruct for models
iter = 10;  % To increase trials, repeat the human trials X times
dstruct = getEmptyDstruct_realData(dstruct_real,iter);
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
        % Re-write choice
        [RDVlist, dstruct(s).choice(ti),dstruct(s).rt(ti),dstruct(s).fixitem{ti},dstruct(s).fixdur{ti},dstruct(s).tItem(ti,:),~] = run_aDDM_rdv(dstruct(s).itemval(ti,:),aGamma_k,dt,k,sig2_k,decbound,yseq_toMaxTime);
        
        if length(dstruct(s).fixdur{ti}) == 6
            figure; hold on;
            x_intervals = [dstruct(s).fixdur{ti}];
            y = [-1, -1, 1, 1]; 
            start_x = 0;                    
            for i = 1:length(x_intervals)                        
                end_x = start_x + x_intervals(i); 
                x = [start_x*1000, end_x*1000, end_x*1000, start_x*1000];                     
                if mod(i, 2) == 1
                    fill(x, y, [0.74, 0.92, 0.99], 'FaceAlpha', 0.3, 'EdgeColor','none');  %Draw blue shading in odd intervals
                else
                    fill(x, y, [1.00,0.77,0.91], 'FaceAlpha', 0.3, 'EdgeColor','none');
                end                        
                start_x = end_x;
            end
            plot(1:length(RDVlist), RDVlist, 'LineWidth', 2, 'Color', [0.00, 0.45, 0.74]);
            xlabel('Time (ms)');
            ylabel('Relative Decision Value');
            title('Choice = ');
            count = count + 1;
            xlim([0, 1000 * sum(x_intervals)]);  
            ylim([-1, 1]); 
            text(1000, 1, sprintf('RV scenario value = %d; FV scenario value = %d', dstruct(s).itemval(ti,1), dstruct(s).itemval(ti,2)));
        end
        if count > 3 %three drawings of the immediate situation
            break;
        end
    end
end

