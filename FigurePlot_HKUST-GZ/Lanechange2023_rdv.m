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

%% Show momentary evidence across time
z = [3,3];  % Item values
sig2_x = 0.00057;  % Evidence accumulation noise
aGamma = 0.3; 
dt = 0.001;  % Time step
t = 3.473;  % Length of time to show
ts = dt*1000:dt*1000:t*1000; 
N = t/dt;
colors_item = [0.34, 0.53, 0.70; 0.91, 0.28, 0.28];

fixseq = [ones(1,round(0.288*N)),2*ones(1,round(0.144*N)),ones(1,round(0.1728*N)),2*ones(1,round(0.222*N)),ones(1,round(0.066*N))];
fixseq = [fixseq,2*ones(1,N-length(fixseq))];
switchpts = ts([0,diff(fixseq)~=0]==1);

% Plot evidence distribution across time 
dx = nan(N,2);   % Evidence accumulation
tItem = [0,0];   % Time points each item was attended to
for n = 1:N
    y = fixseq(n);  % Currently attended item
    tItem(y) = tItem(y)+dt;  % Time advances
    % Evidence accumulation
    dx1 = randn*sqrt( (sig2_x*(aGamma^(2-y)))*dt ) + z(1)*dt; %FV
    dx2 = randn*sqrt( (sig2_x*(aGamma^(y-1)))*dt ) + z(2)*dt; %RV
    dx(n,:) = [dx1,dx2];
end
% Plot
figure; hold on;
for i = 1:2
    if i==1, plots(i) = plot(ts,dx(:,i),'o','Color',colors_item(i,:),'markersize',3);
    else, plots(i) = plot(ts,dx(:,i),'.','Color',colors_item(i,:),'markersize',6);
    end
    plot(get(gca,'xlim'),[z(i)*dt,z(i)*dt],'--','Color',colors_item(i,:));
end
set(gca,'ylim',[-0.002,0.008]);
for sw = 1:length(switchpts), plot([switchpts(sw),switchpts(sw)],get(gca,'ylim'),'k--'); end
xlabel('Time(ms)');
ylabel('Momentary evidence');
legend(plots,{'FV','RV'});
