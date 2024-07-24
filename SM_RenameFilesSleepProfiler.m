function SM_RenameFilesSleepProfiler(source_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT EXAMPLE:
% source_path = fullfile('D:\Studies\05_SleepProfiler\01_SourceData');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(source_path)

% Define the BIDS IDs, subject IDs, age, sex, and education arrays
bids_id = ['sub-SP001'; 'sub-SP002'; 'sub-SP003'; 'sub-SP004'; 'sub-SP005'; ...
     'sub-SP006';'sub-SP007'; 'sub-SP008'; 'sub-SP009'; 'sub-SP010'; ...
     'sub-SP011'; 'sub-SP012'; 'sub-SP013'; 'sub-SP014'; 'sub-SP015'; ...
     'sub-SP016'; 'sub-SP017'; 'sub-SP018'; 'sub-SP019'; 'sub-SP020'; ...
     'sub-SP021';'sub-SP022'; 'sub-SP023'; 'sub-SP024'; 'sub-SP025'; ...
     'sub-SP026'; 'sub-SP027'; 'sub-SP028'; 'sub-SP029'; 'sub-SP030'; ...
     'sub-SP031'; 'sub-SP032'; 'sub-SP033'; 'sub-SP034'];

id = [{'AA_3R079_L'}; {'AA_3R152_L'}; {'AA_3R424_L'};{'AA_3R438_L'}; ...
    {'AA_3R457_L'}; {'AA_3R487_L'}; {'AA_4R143_L'};{'AA_4R342_L'}; ...
    {'COV219'};  {'COV230'}; {'COV250'};  {'COV251'}; {'COV263'};...
    {'COV282'}; {'COV288'};{'COV290'}; {'COV291'} ; {'COV292'};...
    {'COV294'} ;{'COV298'}; {'COV302'};  {'COV309'}; {'COV312'};...
    {'COV317'};  {'COV321'};  {'COV327'}; {'COV_R006'}; {'COV_R027'}; 
    {'COV_R028'}; {'COV_R029'}; {'COV_R043'}; {'COV_R046'};...
    {'COV_R074'}; {'COV_R075'}];


% Get list of folders in source path
Folders = dir(source_path);

% Loop through each folder and rename files
for ii = 1:numel(id)
    % Get subject's old ID and BIDS ID
    sub_OldId = char(id(ii,:));
    sub_BidsId = bids_id(ii, :); 

    % Get the list of files that match the old ID
    OldIdFiles = dir(fullfile(source_path, '**', [sub_OldId '*'])); 

    % Rename the files to use the BIDS ID
    renameFiles(OldIdFiles, sub_BidsId, sub_OldId);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to rename files
function renameFiles(OldIdFiles, newId, oldId)
%EXAMPLE inputs:
% OldIdFiles = OldIdFiles;
% newId = sub_BidsId;
% oldId = sub_OldId;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:numel(OldIdFiles)
        oldFileName = fullfile(OldIdFiles(jj).folder, OldIdFiles(jj).name);
        [pathstr, name, ext] = fileparts(oldFileName);
        newFileName = fullfile(pathstr, strrep(name, oldId, newId));
        movefile(oldFileName, [newFileName, ext]);
    end
end