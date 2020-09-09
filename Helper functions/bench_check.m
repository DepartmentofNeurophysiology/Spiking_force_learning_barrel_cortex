function [ parameters_out] = bench_check( param, iter, scale )
% Performe a benchmark test on the validate set of trials. Do this by
% defining random output weights BPhi, scaled to the learned weights BPhi.
%% Set parameters
G = param.par_comb(iter,1); Q = param.par_comb(iter, 2); 
N = param.N;  Ein = param.Ein; rate = param.rate;
N_train = param.N_train; N_total = param.N_total; N_test = param.N_test;
k = 1; % number of outputs
%% Initialise weights: checken of dit wel dezelfde iedere keer
rng(0) 
p = 0.1; %Set the network sparsity 
E = (2*rand(N,k)-1)*Q;  % Feedback weights
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p); % Static weights
for i = 1:1:N
    QS = find(abs(OMEGA(i,:))>0);
    OMEGA(i,QS) = OMEGA(i,QS) - sum(OMEGA(i,QS))/length(QS);
end
%% Define some random output weights, scale them with the input param.
BPhi = 0.1*(rand(N,k) - 1); 
%% Save weights to structure which will be input to the network
weights.output_weights = BPhi;
weights.static_weights = OMEGA;
weights.feedback_weights = E;
%% Get all the file names of the spiking structures
filePattern = fullfile('..\Spiking structures', '*.mat'); 
theFiles = dir(filePattern);
names = {theFiles.name}; 
shuffled_names = names(randperm(length(names))) ; % Shuffle the file names
%% Make train set and validation set
train = 1:1:N_train * N_total; % Trials to train on.
val = length(train) + 1:1:length(train) +1 + N_test; % Trials to validate on.
%% Benchmark part
save_val_trials = {};
FORCE = false;
for trial = 1: N_test % parfor loop hier?
    disp(['Validating trial = ', num2str(trial)])
    trial_val = shuffled_names(val(trial)); % Select from the test set.
    save_val_trials{trial} = trial_val;
    load( ['..\Spiking structures\',trial_val{1,1}]) % Load the spiking structure from the spiking directory
    pole = SpikeTrainStruct{1,1}.ytrain;
    % Get spiking input and target funtion for a specfic trial
    [ input_struct ] = reservoir_input( SpikeTrainStruct, 1, Ein, N, pole, rate );
    %% MAIN NETWORK FUNCTION   
    [val_dat{trial,1}, val_dat{trial,2}, val_dat{trial,3}, val_dat{trial,4}, val_dat{trial,5}] = spiking_network( param, weights, input_struct, FORCE, iter);
    
    [stats{trial,1}, stats{trial,2}, stats{trial,3}] = spike_stats( val_dat{trial,5}, length(val_dat{trial,4}) , N );
end
%% Save data    
parameters_out.error = val_dat(:,1);
parameters_out.BPhi = val_dat(:,2);
parameters_out.target = val_dat(:,3);
parameters_out.network = val_dat(:,4);
parameters_out.spikes = val_dat(:,5);
parameters_out.A_t= stats(:,1);
parameters_out.ISI= stats(:,2);
parameters_out.Cv= stats(:,3);
parameters_out.val_trials = save_val_trials;
[ acc ] = val_acc( parameters_out);
parameters_out.val_acc = acc;
end

