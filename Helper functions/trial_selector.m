function [shuffled_train, shuffled_val] =...
    trial_selector(prox_touch, dist_no_touch, N_train, N_test)
% trial_selector gets the shuffled and random train and test trials 
% in the ratio 1:1, prox:dist
% Input:
%   * prox_touch = struct containing data regarding proximal whisker touches
%   * dist_no_touch = struct containing data regarding distal whisker touches
%   * N_train = number of training trials
%   * N_test = number of test trials
% Output:
%   * shuffled_train = struct containing shuffled training trials
%   * Shuffled_val = struct containing shuffled test trials

rng('shuffle')

% create indices lists for both conditions
prox_index = 1 : length(prox_touch);
dist_index = 1 : length(dist_no_touch);

% shuffle the proximal and distal trials
prox_index = prox_index(randperm(length(prox_index))); 
dist_index = dist_index(randperm(length(dist_index))); 

% assure the train trials have the correct ratio (1:1, prox:distal)
train_prox = prox_index(1: N_train/2); 
train_dist = dist_index(1: N_train/2);

% collect the data from the structs
index = 1;
for i = 1 : length(train_prox)
    
    % get the struct field
    train_all(index) = prox_touch(train_prox(i));
    
    index = index + 1;
end

for i = 1 : length(train_dist)
    
    % get the struct field
    train_all(index) = dist_no_touch(train_prox(i));
    
    index = index + 1;
end
    

% select validation trials outside of train trials in correct ratio
val_prox = prox_index( N_train/2 + 1 : N_train/2 + N_test/2);
val_dist = dist_index( N_train/2 + 1 : N_train/2 + N_test/2);

% collect the data from the structs
index = 1;
for i = 1 : length(val_prox)
    
    % get the struct field
    val_all(index) = prox_touch(val_prox(i));
    
    index = index + 1;
end

for i = 1 : length(val_dist)
    
    % get the struct field
    val_all(index) = dist_no_touch(val_prox(i));
    
    index = index + 1;
end

% convert to cells and shuffle
train_cell = struct2cell(train_all); 
val_cell = struct2cell(val_all); 

% shuffle the cells
train_cell = train_cell(:, :, randperm(length(train_all)));
val_cell = val_cell(:, :, randperm(length(val_all))); 

% convert them back to structs
fields = ["trial", "session", "spike_struct", "ytrain", "first_touch", "pole_times"];
shuffled_train =  cell2struct(train_cell, fields, 1);
shuffled_val = cell2struct(val_cell, fields, 1);
end

