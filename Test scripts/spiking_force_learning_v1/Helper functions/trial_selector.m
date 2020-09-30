function [shuffled_train, shuffled_val] =...
    trial_selector(prox_touch, dist_no_touch, N_train, N_test)
%TRIAL_SHUFFLE Summary of this function goes here
%   Detailed explanation goes here

% check if the max amount of trials is exceeded?

% shuffle the proximal and distal trials
prox = prox_touch(randperm(length(prox_touch))); 
dist = dist_no_touch(randperm(length(dist_no_touch))); 

% assure the train trials have the correct ratio (1:1, prox:distal)
train_prox = prox(1:N_train/2); 
train_dist = dist(1:N_train/2);
train_all = [train_prox train_dist];

% select validation trials outside of train trials in correct ratio
val_prox = prox( N_train/2 + 1 : N_train/2 + N_test/2);
val_dist = dist( N_train/2 + 1 : N_train/2 + N_test/2);
val = [val_prox val_dist];

% shuffle the trials
shuffled_val = val(randperm(length(val))) ; 
shuffled_train = train_all(randperm(length(train_all))) ; 

end

