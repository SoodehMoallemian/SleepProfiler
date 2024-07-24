datafunction SM_main_eeg_preprocess(sub_id,session_number, source_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function prepares EEG raw data, performs preprocessing, and saves
% files needed for further investigation (e.g., sleep spindles detection).
% To run this function, eeglab must be installed. The function will save
% the generated files in the subject's folder. 
% You also have to add the dependencies to the path of your MATLAB 
% The script works only for the SLEEPPROFILER headband.
% INPUTS:
%   sub_id: Subject ID
%   source_path: Path to the source folder containing EEG data
%   dep_path: path to the dependencies
%
% EXAMPLE INPUTS:
%   sub_id = 'sub-SP001';
%   session_number = 1;
%   source_path = 'D:\Studies\05_SleepProfiler\02_RawData\Sleep_Profiler\';
%
% You can modify this function's paths for your installation.
%
% by Soodeh Moallemian, Ph.D., Rutgers University-2024.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP-01: Convert EDF to SET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-01: Convert EDF to SET')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
sub_session = sprintf('ses-0%d',session_number);
sub_fold = fullfile(source_path, sub_id, sub_session);

cd(sub_fold)
eeglab

%read the edf file and load it to eeglab
edf_file = spm_select('FPList',sub_fold, '.*\.edf');
edf_file = pop_biosig(edf_file);

% get the frequency rate for future calculations
Fs = edf_file.srate;

%read epochs from the txt report file in a table
hypno_file = spm_select('FPList', fullfile(sub_fold), '.*\.csv');
hypno_tab = readtable(hypno_file);

% Read the duration of each epoch which is 30 seconds
% Define the epochs based on the sleep stages in txt report
%it is 30 seconds

% epoch_duration = 30*ones(size(hypno_tab,1),1); 

%  0 – wake
%  1 – N1
%  2 – N2
%  3 – N3
%  4 – L2 (light N2)
%  5 – REM
%  6 – NOS (Sleep not otherwise specified)
%  9 – INVALID
sleep_stages = hypno_tab.PrimaryAutoStage;
% Find the first and last epoches
first_epoch = find(sleep_stages == 1, 1, 'first');
last_epoch = find(sleep_stages == 1, 1, 'last');
sleep_stages = sleep_stages(first_epoch:last_epoch);

% Exclude epochs that are marked as wake (0) or invalid (9)
valid_sleep_stages = sleep_stages(sleep_stages ~= 0 & sleep_stages ~= 9);
% Create the epoch duration vector with 30 in each element
epoch_duration = 30 * ones(size(valid_sleep_stages, 1), 1);
% Define a name for the raw file and create a .set file (readable for eeglab)
filename = [sub_id '_' sub_session '_RawSleep.set'];
% Define Channel names for DREEM3
EEG = pop_select(edf_file, 'channel', {'LEOG', 'REOG', 'E3', 'HMov', 'HPos', 'Pulse', 'Snore'});
% Save the updated file
EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             STEP-02: Clip EDF based on sleep stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-02: Clip EDF based on sleep stages')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
% First data is calculated based on the frequency rate (FS)
% * the duration of each epoch(which is 30) using the real first epoch 
% (before adding the 1)
first_data = (first_epoch) * Fs * epoch_duration(1);
%last data is calculated based on the frequency rate (FS)
% * 30(which is the duration of each epoch) using the last epoch 
last_data = last_epoch * Fs * epoch_duration(1);
EEG = pop_select(EEG, 'point', [first_data last_data]);
EEG.sleepstage_labels = valid_sleep_stages;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-03: Map sleep stages to all samples of sleep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-03: Map sleep stages to all samples of sleep')
disp('-------------------------------------------------')
disp('-------------------------------------------------')

sleep_stages2 = repmat(valid_sleep_stages, 1, Fs * 30);
EEG.stagenames = reshape(sleep_stages2',size(sleep_stages2,1)*size(sleep_stages2,2),1);
stagenames = containers.Map({0, 1, 2, 3, 4, 5, 6},{'W', 'N1', 'N2', 'N3','L2', 'REM', 'NOS'});
temp = cellfun(@(x) stagenames(x), num2cell(EEG.stagenames), 'UniformOutput', false);
EEG.stagenames = temp;

filename = sprintf('%s_%s_SnippedSleep_filtered_bandpass.set', sub_id, sub_session);
EEG = pop_saveset(EEG, 'filename', filename, 'filepath', sub_fold, 'check', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-04: Visual quality control of the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-03: Visual quality control of the data')
disp('Here you can check the quality of the signals from each channel.')
disp('you can trim(clip) the data vertically or horizentally.')
disp(['HORIZENTALLY: Please visually check the data and if you find one' ...
    'channel that looks off in more than 50 %% of the cases, ' ...
    'delete that channel by simply clining on the name of the channel ' ...
    'on the left side of the figure. '])
disp(['VERTICALLY: When you visually check the signals and ' ...
'using "1" (the start point of the trim) and "2" the end point'])

disp('-------------------------------------------------')
disp('-------------------------------------------------')
EEG = pop_loadset(fullfile(sub_fold,[sub_id '_' sub_session '_SnippedSleep_filtered_bandpass.set']));
load(fullfile (dep_path , 'Dependencies', 'chanlocs_dreem3_mock.mat'), 'chanlocs1');
EEG.chanlocs = chanlocs1;
EEG.original_chanlocs = EEG.chanlocs;
EEG = csc_eeg_plotter_NEW_Eisai(EEG);
EEG.bad_channels{1} = EEG.hidden_channels;
EEG = pop_select(EEG, 'nochannel', EEG.bad_channels{1});
% save original bc data
EEG.save_dir= sub_fold;
EEG.filename = strrep(EEG.filename,'.set','_bc1pass.set');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 STEP-04: Create N2N3 set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-04: Create N2N3 set')
disp('-------------------------------------------------')
disp('-------------------------------------------------')

N2N3_samples_inds = [];
N2only_inds=[];
N3only_inds=[];

for i =1:length(EEG.stagenames)
    if strcmp(EEG.stagenames{i},'N2')||strcmp(EEG.stagenames{i},'N3')
        N2N3_samples_inds(i)=1;
    end
    if strcmp(EEG.stagenames{i},'N2')
        N2only_inds(i)=1;
    end
    if strcmp(EEG.stagenames{i},'N3')
        N3only_inds(i)=1;
    end
end
t_n2n3=find(N2N3_samples_inds); % find the indices both for N2 and N3
dif=find(t_n2n3(2:end)-t_n2n3(1:end-1)>1);
ends=[t_n2n3(dif),t_n2n3(end)];
starts=[t_n2n3(1),t_n2n3(dif+1)];
n2n3_intervals=[starts',ends'];

N2N3_EEG=pop_select( EEG, 'point',n2n3_intervals );

newfilename  = [strrep(N2N3_EEG.filename,'_bc1pass.set','_n2n3_noarousals'),'.set'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       STEP-05: Bad segment detection (semi-automatic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-05: Bad segment detection (semi-automatic)')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('********')
disp('This part of the function is semi-automatic.')
disp('********')
fprintf('Please choose the method you want to use for segmenting the bad signal:\n')
method = 'wispic';
N2N3_EEG=csc_artifact_rejection_automated_Eisai(N2N3_EEG, method, 'epoch_length', 6);
N2N3_EEG = pop_select(N2N3_EEG, 'notime', N2N3_EEG.bad_regions);
N2N3_stages=N2N3_EEG.stagenames(find(N2N3_samples_inds));
N2N3_artifactsamps = zeros(1,length(N2N3_stages));
bad_regions = N2N3_EEG.bad_regions;

for i = 1:length(bad_regions)
    start = bad_regions(i,1)*N2N3_EEG.srate;
    last = bad_regions(i,2)*N2N3_EEG.srate;
    if start == 0
        start = 1;
    end
    N2N3_artifactsamps(start:last) = 1;
end

N2N3_stages_afterfft = N2N3_stages(N2N3_artifactsamps ==0);
cd(N2N3_EEG.save_dir);
save(sprintf('%s_N2N3_stagenames_afterfft.mat',sub_id), ...
    'N2N3_stages_afterfft');
%save file with fft
tic
fftfilename = [strrep(newfilename, '_n2n3_noarousals.set', ...
    '_n2n3_noarousals_fftstep'),'.set'];
EEG = pop_saveset(N2N3_EEG,'filename', fftfilename, 'filepath', ...
    sub_fold, 'check', 'on');
disp(['... elapsed time ',num2str(toc/60),' saving set file']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-06: Check data after artifact rejection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-06: Check data after artifact rejection')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('this is your last chance to look into the data and check its quality visually.')
disp(['Note! Now, on the figure, you will only see the N2N3 sleep stages ' ...
    'which represent the nonREM sleep time course. Therefore, it is normal to see ' ...
    'less on the horizental timeline.'])
disp(['if you have not done any triming before, it is best to add two clips' ...
    ' at the end of the signals (for horizental clipping) '])
N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
if ~isempty(N2N3_EEG.hidden_channels)
    N2N3_EEG.bad_channels{1} = cat(2, N2N3_EEG.bad_channels{1}, ...
        str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels}));
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
    disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      STEP-07: Check the Clippings(they must be in even numbers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-07: Check the Clippings' )
disp('-------------------------------------------------')
disp('-------------------------------------------------')
clip_check=0;
while clip_check == 0
    if isfield(N2N3_EEG,'csc_event_data')
        if ~isempty(N2N3_EEG.csc_event_data)
            if size(N2N3_EEG.csc_event_data,1) > 2
                for n = 2:size(N2N3_EEG.csc_event_data,1)
                    if isequal(N2N3_EEG.csc_event_data(n,3), N2N3_EEG.csc_event_data(n-1,3))
                        disp('there are back to back same color markings');
                        beep
                        N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
                        clip_check=0;
                    else
                        clip_check=1;
                    end
                end
            end
        end
    else
    end
    check_csc_event_data =  ~mod(size(N2N3_EEG.csc_event_data,1),2);

    if check_csc_event_data == 0
        disp('there is an uneven number of markings');
        beep
        N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
        clip_check=0;
    else
        clip_check=1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      STEP-08: Remove epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-08: Remove epochs' )
disp('-------------------------------------------------')
disp('-------------------------------------------------')

if clip_check == 1 && ~isempty(N2N3_EEG.csc_event_data)
    event_starts = cellfun(@(x) strcmp(x, 'event 1'), N2N3_EEG.csc_event_data(:, 1));
    
    % sanity check for artifact event markers
    if sum(event_starts) ~= sum(~event_starts)
        fprintf('\nWarning: uneven number of events, check event_data\n');
    end
    
    % use EEGLAB to remove the points
    N2N3_EEG.bad_segments{1} = [cell2mat(N2N3_EEG.csc_event_data(event_starts, 2)), ...
        cell2mat(N2N3_EEG.csc_event_data(~event_starts, 2))];
    
    % convert the timing from seconds to samples
    N2N3_EEG.bad_segments{1} = floor(N2N3_EEG.bad_segments{1} * N2N3_EEG.srate);
    
    % use EEGLAB to remove the regions
    N2N3_EEG = pop_select(N2N3_EEG, 'nopoint', N2N3_EEG.bad_segments{1});
    
    
end
% save the file with channels/epochs deleted & marked in n2n3 file
tic
n2n3filename = [strrep(fftfilename,'_n2n3_noarousals_fftstep.set','_n2n3_noarousals_fftstep_artifcorr'),'.set'];
N2N3_EEG = pop_saveset(N2N3_EEG,'filename', n2n3filename, 'filepath', N2N3_EEG.save_dir, 'check', 'on');

load(sprintf('%s_N2N3_stagenames_afterfft',sub_id));
numsamps_afterfft = zeros(1,length(N2N3_stages_afterfft));
bad_segments = cell2mat(N2N3_EEG.bad_segments);

for j = 1:size(bad_segments,1)
    numsamps_afterfft(bad_segments(j,1):bad_segments(j,2)) = 1;
end

N2N3_stages_afterfft_snipped = N2N3_stages_afterfft(numsamps_afterfft ==0);
% Removing the last sample, since this is a mismatch, can recheck this
% later
if size(N2N3_stages_afterfft_snipped,2) ~= size(N2N3_EEG.data,2)
    N2N3_stages_afterfft_snipped(end) = [];
end
N2N3_EEG.N2N3_stagenames_afterfft_snipped =  N2N3_stages_afterfft_snipped;
save([sub_fold, '/', sub_id '_N2N3_stagenames_afterfft_snipped.mat'],'N2N3_stages_afterfft_snipped');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-09: Check bad channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-09: Check bad channels')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
finished = 'N';
while ~isempty(N2N3_EEG.hidden_channels) || strcmp(finished,'N')
    N2N3_EEG = csc_eeg_plotter_NEW_Eisai(N2N3_EEG);
       
    new_bad_channels = str2double({N2N3_EEG.chanlocs(N2N3_EEG.hidden_channels).labels});
    
    tic
    N2N3_EEG.bad_channels{1} = cat(2,N2N3_EEG.bad_channels{1},new_bad_channels);
    if ~isempty(N2N3_EEG.hidden_channels)
        N2N3_EEG = pop_select(N2N3_EEG, 'nochannel', N2N3_EEG.hidden_channels);
        disp(['... elapsed time ',num2str(toc/60),' removing bad channels']);
%       N2N3_EEG = pop_reref( N2N3_EEG, [] );
        N2N3_EEG.hidden_channels = [];
    end
   if isfield(N2N3_EEG,'psd'); N2N3_EEG = rmfield(N2N3_EEG, 'psd');end
    for n = 1:N2N3_EEG.trials
        [N2N3_EEG.psd.data(:,:,:,n),N2N3_EEG.psd.Hzbins] = psddata(squeeze(N2N3_EEG.data(:,:,n)),N2N3_EEG.srate,6,40,1); % only use fast for massive data
      
    end
    close all;
    plot_bands_spectra_EEG_PAMMS_Eisai(N2N3_EEG);
    corr_ans = 0;
    while ~corr_ans
        finished = input('Finished Y/N? [Y]','s');
        if ~ismember(finished, {'Y','N'})
            corr_ans = 0;
        else
            corr_ans = 1;
        end
    end
end

%save file data set
tic
if finished == 'Y'
    finishedfilename = [strrep(n2n3filename,'_n2n3_noarousals_fftstep_artifcorr.set','_n2n3_noarousals_fftstep_artifcorr_final'),'.set'];
    %           finishedfilename = [strrep(EEG.filename,'_n2n3_noarousals_artifcorr_fftstep.set','_n2n3_noarousals_artifcorr_fftstep_final'),'.set'];
    N2N3_EEG.datfile = finishedfilename;
    N2N3_EEG = pop_saveset(N2N3_EEG,'filename', finishedfilename, 'filepath', sub_fold, 'check', 'on');
end
disp(['... elapsed time ',num2str(toc/60),' saving set file']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          STEP-10: Segmenting the N2 and N3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beep
disp('-------------------------------------------------')
disp('-------------------------------------------------')
disp('STEP-10: Segmenting the N2 and N3')
disp('-------------------------------------------------')
disp('-------------------------------------------------')
% We will only save the final
% data after checking the atrifact rejection from N2 and N3 seperately

load(fullfile(sub_fold,[sub_id '_N2N3_stagenames_afterfft_snipped.mat'])); 
N2_EEG = N2N3_EEG;

finishedfilename = 'preprocessed_eeg';
N2cleaned_filename = [sub_id '_' finishedfilename '_N2only.set'];
N2only_cleaned_ind = (strcmp(N2N3_stages_afterfft_snipped,'N2'));

swa_selectStagesEEGLAB_BAR(N2_EEG, N2only_cleaned_ind,'N2only_cleaned',N2cleaned_filename);
% N2only_cleaned = pop_saveset(N2_EEG, 'filename', N2cleaned_filename, 'filepath', sub_fold ,'check', 'on');
 

N3_EEG = N2N3_EEG;
N3cleaned_filename = [sub_id '_' finishedfilename '_N3only.set'];

N3only_cleaned_ind = (strcmp(N2N3_stages_afterfft_snipped,'N3'));
swa_selectStagesEEGLAB_BAR(N3_EEG, N3only_cleaned_ind,'N3only_cleaned',N3cleaned_filename);
% N3only_cleaned = pop_saveset(N3_EEG, 'filename', N3cleaned_filename, 'filepath',sub_fold ,'check', 'on');



NREMpsd = N2N3_EEG.psd;
save NREMpsd NREMpsd;

end