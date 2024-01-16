clear;

rootdir = 'C:/Users/A/Documents/MATLAB/driving_model/ActiveSensing';
addpath(genpath(rootdir));
data_dir = fullfile(rootdir,'data');
set(0,'defaultAxesFontSize',10);

load(fullfile(data_dir,sprintf('datastruct_hkust0902.mat')));
dstruct = dstruct_real;

%% Show fixation bias effects
showfixationdurations = 1;  
showBehavPlots = 1; 
showFixBiasEffects = 1; 

nbin_fn = 4;
nbin_valsum = 5;
useFixationBins = 1;

% Show fixation bias effects
if showFixBiasEffects==1
    nbin_fn = 4;
    nbin_valsum = 5;
    useFixationBins = 1;
end

% Show distribution of fixation durations
if showfixationdurations==1
    numfixshow = 8; 
    % fixdur_all = cell(1,numfixshow);
    for s = 1:length(dstruct)
        for fn = 1:numfixshow
            fixdur_sub_all = [];
            for ti = 1:length(dstruct(s).fixdur)
                if length(dstruct(s).fixdur{ti}) >= fn
                    fixdur_sub_all = cat(1,fixdur_sub_all,dstruct(s).fixdur{ti}(fn));
                end
            end
            fixdur_sub(s,fn) = mean(fixdur_sub_all);
        end
    end

    [mean_fixdur,se_fixdur] = getMeanAndSE(fixdur_sub);

    % Add boxplot
    figure; hold on;
    bp = boxchart(fixdur_sub);
    set(bp, 'BoxEdgeColor', 'black'); 
    set(bp,'BoxFaceColor', [0.529, 0.808, 0.922] ); 
    % Add mean
    mean_values = mean(fixdur_sub);
    scatter(1:numfixshow, mean_values, 'x', 'k', 'SizeData', 100, 'LineWidth', 1);  % Add "x" to the mean position and adjust the size and thickness
    ylabel('Fixation duration (s)'); xlabel('Fixation number'); pbaspect([1,1,1]);
    xticklabels(1:numfixshow); 
end


% Show distribution of fixation durations for choice=1
if showfixationdurations==1
    numfixshow = 7; %
    % fixdur_sub_choice1 = zeros(length(dstruct), numfixshow);
    for s = 1:length(dstruct)
        for fn = 1:numfixshow
            for row = 1:length(dstruct(s).choice)
                if dstruct(s).choice(row,1) == 1
                    fixdur_sub_all_choice1 = [];
                    %for ti = 1:length(dstruct(s).fixdur)
                    if mod(fn,2) == 0
                        if length(dstruct(s).fixdur{row}) >= fn
                             fixdur_sub_all_choice1 = cat(1,fixdur_sub_all_choice1,dstruct(s).fixdur{row}(fn));
                        end
                        %end
                        fixdur_sub_choice1(s,fn) = mean(fixdur_sub_all_choice1);
                    end
                end

            end
        end
    end

    fixdur_sub_choice1 = fixdur_sub_choice1(:, mod(1:size(fixdur_sub_choice1, 2), 2) == 0);
    [mean_fixdur_choice1,se_fixdur_choice1] = getMeanAndSE(fixdur_sub_choice1);

% Show distribution of fixation durations for choice=2

    for s = 1:length(dstruct)
        for fn = 1:numfixshow
            for row = 1:length(dstruct(s).choice)
                if dstruct(s).choice(row,1) == 2
                    fixdur_sub_all_choice2 = [];
                    %for ti = 1:length(dstruct(s).fixdur)
                    if mod(fn,2) == 0
                        if length(dstruct(s).fixdur{row}) >= fn
                             fixdur_sub_all_choice2 = cat(1,fixdur_sub_all_choice2,dstruct(s).fixdur{row}(fn));
                        end
                        %end
                        fixdur_sub_choice2(s,fn) = mean(fixdur_sub_all_choice2);
                    end
                end

            end
        end
    end

    fixdur_sub_choice2 = fixdur_sub_choice2(:, mod(1:size(fixdur_sub_choice2, 2), 2) == 0);
    [mean_fixdur_choice2,se_fixdur_choice2] = getMeanAndSE(fixdur_sub_choice2);

    figure; hold on;
    Choise = repelem([1; 2], 12);%Repeat the matrix 12 times in the vertical direction
    Fixnum = repelem([1; 2; 3], 4);
    Fixnum = repmat(Fixnum, 2, 1);
    fixdur_sub_choice = cat(1, fixdur_sub_choice1.', fixdur_sub_choice2.');
    fixdur_sub_choice = reshape(fixdur_sub_choice',[],1);

    bc_all = boxchart(Fixnum, fixdur_sub_choice, 'GroupByColor', Choise);
    set(bc_all, 'BoxEdgeColor', 'black'); 
    ylabel('Fixation duration (s)');xlabel('Fixation number');
    legend('Lane change' , 'Go straight');
end