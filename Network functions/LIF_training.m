function [ training_output ] = LIF_training_v1(param, scale_param, savefolder)
%LIF_TRAINING Summary of this function goes here
%   Detailed explanation goes here



%% Parameters

% input parameters
N = param.N;
N_th = param.N_th;
rate = param.rate;
N_train = param.N_train;
N_total = param.N_total;
N_test = param.N_test;

% weight scaling parameters
Win = scale_param.Win;
G = scale_param.G;
Q = scale_param.Q;

% set the static network sparsity
p = 0.1;

% set input network sparsity
Winp = scale_param.Winp;

% number of network outputs 
k = 1;

%% Initialize matrix weights

% input weights
input = Win.*rand(N,N_th).*(rand(N,N_th) < Winp); 

rng('shuffle')

% static weights

static =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);

% set the row average to be zero, explicitly 
for i = 1:1:N 
    QS = find(abs(static(i,:))>0);
    static(i,QS) = static(i,QS) - sum(static(i,QS))/length(QS);
end

%{
% Apply Dale's law

% number of exci and ini neurons 
Nexc = round(0.5 * N);
Nini = N - Nexc;

% initiate both halves of the static weights
static_exc = G*abs(randn(N, Nexc)).*(rand(N, Nexc) < p)/(sqrt(N)*p);
static_ini = G*abs(randn(N, Nini)).*-(rand(N, Nini) < p)/(sqrt(N)*p);

% set the mean input of each neuron to zero by adjusting the inibition
for i = 1:N 
    QS = find(static_ini(i, :) ~= 0);
    
    row_sum = sum(static_exc(i, :)) + sum(static_ini(i, :));
    static_ini(i, QS) = static_ini(i, QS) - row_sum/length(QS);
end

static = [static_exc static_ini];
%}

% feedback weights
feedback = Q*(2*rand(N,k)-1);

% output weights 
output = zeros(N, k);

% store all weights
weights.static = static;
weights.input = input;
weights.feedback = feedback;
weights.output = output;

%% list of trials for training and validation
train_trials = param.train_trials;
test_trials = param.test_trials;

%% Train and test the network

% run for N amount of epochs
for epoch = 1:N_total
    
    disp(['Epoch nr. ', num2str(epoch)])
    
    
    % TRAINING
    disp(['Training network, number of trials = ', num2str(N_train)])
    
    % apply FORCE learning
    FORCE = param.FORCE;
    
    % weight change storage variable
    weight_change = zeros(N_train, 1);
    
    % training trials
    for trial = 1:N_train
        
        % get the trial name 
        trialId = train_trials(trial).trial;
        
        % make or load the spikes
        if param.makespikes
            
            % get the trial session and create the spikingstruct
            session = train_trials(trial).session;
            SpikeTrainStruct = make_trial_spikes(session, trialId);
        else
            % get the struct name and load it
            trial_mat = train_trials(trial).spike_struct;
            load( ['./Spiking structures/', trial_mat]);
        end
        
        % get the pole location and the input struct and target function
        pole = train_trials(trial).ytrain;
        [thalamus_input, target] =...
            reservoir_input(SpikeTrainStruct, 1, input, N, pole, rate); 
        
        % SIMULATE NETWORK
        % save the old output weights
        old_output = output;
        [~, output, ~, ~, ~] =...
            LIF_spiking_network_v1(param, weights, thalamus_input, target, FORCE);
        
        % calculate the weight difference and update the output weights
        d_output = old_output - output;
        weight_change(trial, 1) = sum(abs(d_output));
        weights.output = output;
    end
    
    
    % VALIDATION
    disp(['Testing network, number of trials = ', num2str(N_test)])
    
    % turn off FORCE learning
    FORCE = false;
    
    % validation trials
    for trial = 1:N_test
        
        % get the trial name 
        trialId = test_trials(trial).trial;
        
        % make or load the spikes
        if param.makespikes
            
            % get the trial session and create the spikingstruct
            session = test_trials(trial).session;
            SpikeTrainStruct = make_trial_spikes(session, trialId);
        else
            % get the struct name and load it
            trial_mat = test_trials(trial).spike_struct;
            load( ['./Spiking structures/', trial_mat]);
        end
        
        % save the validation trials and firs touches
        test_output.trials{trial} = trialId;
        test_output.first_touches(trial,1) = test_trials(trial).first_touch;
        
        % get the pole location and the input struct and target function
        pole = test_trials(trial).ytrain;
        [thalamus_input, target] =...
            reservoir_input(SpikeTrainStruct, 1, input, N, pole, rate);
        
        % SIMULATE NETWORK
        [ error, output_weights, Zx, Z_out, tspikes ] =...
            LIF_spiking_network_v1(param, weights, thalamus_input, target, FORCE);
        
        % get the spikinging statistics
        trial_length = length(Z_out);
        [ A_t, ISI, Cv ] = spike_stats( tspikes, trial_length , N );
        
        weights.output = output_weights;
        test_output.error{trial} = error;
        test_output.Zx{trial} = Zx;
        test_output.Z_out{trial} = Z_out;
        test_output.tspikes{trial} = tspikes;
        test_output.stats{trial}.A_t = A_t;
        test_output.stats{trial}.ISI = ISI;
        test_output.stats{trial}.Cv = Cv;
    end
    
    % calculate the test accuracy
    acc = val_acc(N_test, test_output.Z_out, test_output.Zx, test_output.first_touches);
    
    % save the training data per epoch
    
    training_output(epoch).weight_change = weight_change;
    training_output(epoch).train_trials = train_trials;
    training_output(epoch).test_output = test_output;
    training_output(epoch).acc = acc;
    training_output(epoch).weights = weights;
end

%% save the output of the network

output_save.training_output = training_output;
output_save.scale_param = scale_param;

f = filesep;

filename = '';
field_names = fieldnames(scale_param);
for i = 1 : length(field_names)
    field = char(field_names(i));
    value = getfield(scale_param, field);
    str = [field '_' num2str(value)];
    
    filename = [filename str];
end

 savename = [savefolder f filename '.mat'];
 save(savename, 'training_output', 'scale_param')

