function [L1A,cycles_L1A]=Lib_Load_S3A_L1A_Data_products_Rathbun()



% Define the folder path
folderPath = 'F:\Sentinel-3-Level1\Lakes\317 Rathbun\L1A';

% Get a list of all .nc files in the folder
filePattern = fullfile(folderPath, '**', '*.nc'); % Using '**' for recursive search
files = dir(filePattern);

% Initialize an empty cell array to store the file paths
filePaths = {};
L1A = cell(1, length(filePaths));

% Iterate through the files and store their paths
for k = 1:length(files)
    % Get the full file path
    fullPath = fullfile(files(k).folder, files(k).name);
    
    % Extract the date from the filename (assuming date is in format YYYYMMDD)
    % cycleNumbers = char(regexp(fullPath, '_\d{3}_', 'match'))
        % Create the variable name dynamically
        % varName = ['L1A' cycleNumbers];
            L1A{k} = (fullPath);

end

cycles_L1A = [];
for i = 1:numel(L1A)
    cycleNumbers = char(regexp(string(L1A(i)), '_\d{3}_', 'match'));
    tmp = str2double(cycleNumbers(2:4));
    cycles_L1A = [cycles_L1A; tmp];
end



end