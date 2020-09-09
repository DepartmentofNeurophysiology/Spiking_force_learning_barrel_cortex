function [ theFiles ] = file_names( myFolder )
% Save all filenames in 'myFolder' and return to a struct

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% Base file name is: theFiles(k).name


end

