function [ acc ] = val_acc(test_trials, network_output, target, first_touches)
% Calculate the validation accuracy 

%% Calculate the mean value each trial

mean_network = zeros(test_trials, 1);
mean_target = zeros(test_trials, 1);

for i = 1:test_trials
    z_t = network_output{i};
    zx = target{i};
      
    % trial_length = z_t - 1200;
    % Only calculate mean when output is trained
    %s_t = length(zx) - 800 - 500;
    s_t = first_touches(i);
    if s_t <= 0 
        s_t = 1000;
    end

    %s_t = 1000;
    mean_target(i) =  mean(zx(s_t:end));
    mean_network(i) = mean(z_t(s_t:end));
end
%% Set a classifcation boundary, mean positive value or mean negative value.
network_answer = zeros(test_trials,1);
correct_answer = zeros(test_trials, 1);
for i = 1:test_trials
    if mean_network(i) > 0
        network_answer(i) = 1;
    else
        network_answer(i) = -1;
    end
    if mean_target(i) > 0
        correct_answer(i) = 1;
    else
        correct_answer(i) = -1;
    end
end
%% Calculate the accuracy 
multi = network_answer .* correct_answer; % All positive values indicated correctly classified val_trials
correct = sum(multi==1); % Count total positive values
acc = correct/test_trials;
end

