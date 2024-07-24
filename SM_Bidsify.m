% Define the source directory
sourceFolder = pwd;
% Get a list of all files in the source folder
fileList = dir(fullfile(sourceFolder, 'sub-SP0*_N*'));

% Loop through each file
for k = 1:length(fileList)
    % Get the current file name
    currentFileName = fileList(k).name;
    
    % Extract the subject ID and session number
    tokens = regexp(currentFileName, 'sub-SP0(\d{2}+)_N(\d+)', 'tokens');
    if ~isempty(tokens)
        subID = tokens{1}{1};  % Extracted subject ID
        sessionNumber = tokens{1}{2};  % Extracted session number
        
        % Create the new file name
        newFileName = sprintf('sub-SP0%s_ses-0%s', subID, sessionNumber);
        
        % Create the new folder name
        newFolderName = fullfile(sourceFolder, sprintf('ses-0%s', sessionNumber));
        
        % Create the new folder if it doesn't exist
        if ~exist(newFolderName, 'dir')
            mkdir(newFolderName);
        end
        
        % Move and rename the file
        movefile(fullfile(sourceFolder, currentFileName), fullfile(newFolderName, newFileName));
    end
end

disp('Files have been renamed and moved successfully.');
