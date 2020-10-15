function [ s_correct ] = import_data( myFolder )
% Takes a folder to import data from, and returns the correct data
% seperated into curvature and angle. 

%% Import structures from folder
filePattern = fullfile(myFolder, '*.mat'); 
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  struct = load(fullFileName);
  s(k) = struct;
end
%% Seperate into angle and curve traces and take only correct trials
for k =1:length(theFiles)
    len = length(s(k).xtrain(1,:))/2;
    s(k).x_c = s(k).xtrain(:,1:len) ;
    s(k).x_a = s(k).xtrain(:,len+1:end) ;
    corr_trials = find(s(k).correct == 1);
    s_correct(k).ytrain = s(k).ytrain(corr_trials,:);
    s_correct(k).x_c = s(k).x_c(corr_trials,:);
    s_correct(k).x_a = s(k).x_a(corr_trials,:);
    s_correct(k).trial_number = s(k).trial_number(corr_trials);
    s_correct(k).name = theFiles(k).name;
end


end

