function [ acc ] = val_acc( parameters_out, first_touches )
% Calculate the validation accuracy 
num_val = length(parameters_out.network);
mean_network = zeros(num_val, 1);
mean_target = zeros(num_val, 1);


%% Calculate the mean value each trial
for i = 1:num_val
    zx = parameters_out.network{i};
    z_t = parameters_out.target{i};
    trial_length = z_t - 1200;
    % Only calculate mean when output is trained
    %s_t = length(zx) - 800 - 500;
    s_t = first_touches(i);
    if s_t == 0
        s_t = 1000;
    elseif s_t < 0
        s_t = 1000;
    end
    %s_t = 1000;
    mean_network(i) =  mean(zx(s_t:end));
    mean_target(i) = mean(z_t(s_t:end));
end
%% Set a classifcation boundary, mean positive value or mean negative value.
network_answer = zeros(num_val,1);
correct_answer = zeros(num_val, 1);
for i = 1:num_val
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
acc = correct./num_val;
end

