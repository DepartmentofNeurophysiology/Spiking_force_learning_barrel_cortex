function [shuffled_train, shuffled_val] = left_right_trace(name_list, N_train, N_test, N_total)
%% Left and right
names = name_list;
%
pole = zeros(length(names),1);

for i = 1:length(names)
    pole(i,1) = names{i,4}{1,1}(1);
end

left = find(pole == 1);
right = find(pole == -1);
%

%train_nr = N_train * N_total;
train_nr = N_train;
train_left = left(1:train_nr/2);
train_right = right(1:train_nr/2);

train_all = [train_left; train_right];

val_left = left( train_nr/2 +1 : train_nr/2+ N_test/2);
val_right = right( train_nr/2 +1 : train_nr/2+ N_test/2);

val = [val_left; val_right];

shuffled_val = val(randperm(length(val))) ; % Shuffle 
shuffled_train = train_all(randperm(length(train_all))) ; % Shuffle 

end