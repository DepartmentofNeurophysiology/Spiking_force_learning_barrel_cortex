function [ dat ] = shuffle_data( s_correct )

num_trials = [];

for k = 1:length(s_correct)
    num_trials = num_trials + length(s_correct(k).trial_number);
end

all_data = cell(num_trials, 5);
% Columns: x_c, x_a, ytrain, trial number, struct number
index = 0;
for k = 1:length(s_correct)
    shape = size(s_correct(k).x_c);   
    for i = 1:shape(1)
        all_data{i + index,1} = s_correct(k).x_c(i,:);
        all_data{i + index,2} = s_correct(k).x_a(i,:);
        all_data{i + index,3} = s_correct(k).ytrain(i,:);
        all_data{i + index,4} = s_correct(k).trial_number(i,:);
        all_data{i + index,5} = k;
        all_data{i + index,6} = s_correct(k).name;
    end
    index = index + shape(1);
end

% Shuffle all the rows

shuffled_data = all_data(randperm(size(all_data, 1)), :);

dat.x_c = shuffled_data(:,1);
dat.x_a = shuffled_data(:,2);
dat.ytrain = shuffled_data(:,3);
dat.trial_number = shuffled_data(:,4);
dat.struct_number = shuffled_data(:,5);
dat.struct_name = shuffled_data(:,6);
end

