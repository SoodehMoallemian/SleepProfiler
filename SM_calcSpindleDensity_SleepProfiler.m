function SM_calcSpindleDensity_SleepProfiler (source_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will calculate the mean,duration(seconds) of the spindles 
% within N2 and N3 stages of sleep. It will read the data from hyponogram 
% excel files ...
% The data was acquiried using the SleepProfiler headband.
% INPUT EXAMPLE: source_path = fullfile('D:\Studies\05_SleepProfiler\SleepProfiler_Splitted_Hypnograms')
%
%
% written by Soodeh Moallemian Ph.D. Rutgers university, June 2024 
% s.moallemian@rutgers.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go to the source directory
cd(source_path);

% Find all CSV files ending with 'N1'
files = dir(source_path);

% Initialize an empty table
data_table = [];

% Loop through each file
for i = 3:length(files)
  
  % Read the current CSV file
  data = readtable(files(i).name);
  
  % Select columns AA and BW
  selected_data = data(:, {'PrimaryAutoStage','Sleep_WakeActigraphy',...
      'Spindles','Sigma','SigmaU','SpDurSec','AveAlpMx','AveSigAve',...
      'SpindlesU','SpDurSecU','Spindles1','Spindles16','SpindlesU1',...
      'SpindlesU16'});
  
  % Append the data to the table
  data_table = [data_table; selected_data];
end

% Check if any files were found
if isempty(data_table)
  warning('No CSV files ending with N1 were found in the source directory.');
  return;
end

% Save the table to a new file 
writetable(data_table, [files(i).name(1:end-4) '_SpindleInfo_N1']);

end
