function [fixed_train, fixed_val] =...
    fixed_trial_selector(prox_touch, dist_no_touch, N_train, N_test)
% fixed_trial_selector gets the fixed whisker trials to train and test the
% network, by selecting one proximal and one distal trial 
% Input:
%   * prox_touch = struct containing data regarding proximal whisker touches
%   * dist_no_touch = struct containing data regarding distal whisker touches
%   * N_train = number of training trials
%   * N_test = number of test trials
% Output:
%   * fixed_train = struct containing shuffled training trials
%   * fixed_val = struct containing shuffled test trials

rng(0) % same sequence of random numbers

% randomly selects one distal and one proximal trial 
train_prox(1:N_train/2) = randi([1,length(prox_touch)],1);
train_dist(1:N_train/2) = randi([1,length(dist_no_touch)],1);

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

% selects the same validation trials as for train trials 
val_prox(1:N_test/2) = train_prox(1:N_test/2);
val_dist(1:N_test/2) = train_dist(1:N_test/2);

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
fixed_train =  cell2struct(train_cell, fields, 1);
fixed_val = cell2struct(val_cell, fields, 1);
end

