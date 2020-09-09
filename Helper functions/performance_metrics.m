function [ perform ] = performance_metrics( run, which_run, EPOCHS, VAL_BATCH )
% Returns the errors, the sum of errors for an EPOCH, the sum of weight
% changes per EPOCH, and the mean value of the network during target
% function. 
% Also returns the network answer and the correct answer based on the mean
% values. With this the classification accuracy per epoch can be
% classified.

a = which_run;
errors = [];
sum_error = zeros(EPOCHS,1);
mean_network = zeros(EPOCHS, VAL_BATCH);
mean_target = [];

for i = 1:EPOCHS
    for n = 1:(VAL_BATCH)
        errors = [errors run(a).parameters_out(i).error{n}];
        sum_error(i,1) = sum_error(i,1) + run(a).parameters_out(i).error{n};
        % Network output
        zx = run(a).parameters_out(i).network{n};
        z_t = run(a).parameters_out(i).target{n};
        % Only calculate mean when output is trained
        s_t = length(zx) - 800 - 500;
        mean_network(i,n) =  mean(zx(s_t:end));
        mean_target(i,n) =  mean(z_t(s_t:end));
    end
end

perform.errors = errors;
perform.sum_error = sum_error;
perform.weight_changes = weight_change( run(a).parameters_out);
perform.mean_network = mean_network;
perform.mean_target = mean_target;

network_answer = zeros(EPOCHS, VAL_BATCH);
correct_answer = zeros(EPOCHS, VAL_BATCH);

for i = 1:EPOCHS
    for n = 1:VAL_BATCH
        if mean_network(i,n) > 0
            network_answer(i,n) = 1;
        else
            network_answer(i,n) = -1;
        end
        if mean_target(i,n) > 0
            correct_answer(i,n) = 1;
        else
            correct_answer(i,n) = -1;
        end
    end
end

perform.network_answer = network_answer;
perform.correct_answer = correct_answer;

%% Calculate the accuracy per epoch
multi = network_answer .* correct_answer; % All positive values indicated correctly classified val_trials
correct = sum(multi==1, 2); % Count total positive values
val_acc = correct./VAL_BATCH;
perform.val_acc = val_acc;

end

