function [training_output] = LIF_training(param, scale_param, savefolder)
% LIF_TRAINING prepares the weights and input of the network,
% and runs the test and train trials of the network with spiking neurons,
% based on the leaky integrate and firing model
% Input:
%   * param = input parameters
%   * scale_param = weight scaling parameters
%   * savefolder = string with path to the savefolder
% Output:
%   * training_output = a structure containing the training output

f = filesep;
%% Parameters
% input parameters
N = param.N;                % number of neurons in the reservoir
N_th = param.N_th;          % number of thalamus neurons
N_train = param.N_train;    % number of training trials
N_test = param.N_test;      % number of validation trials
N_total = param.N_total;    % number of epochs
rate = param.rate;          % rate of the intermediate poisson firing (Hz)

% weight scaling parameters
Win = scale_param.Win;      % the input weights
G = scale_param.G;          % the static weights
Q = scale_param.Q;          % the feedback weights
Winp = scale_param.Winp;    % network sparsity
Pexc = scale_param.Pexc;    % percentage of excitatory neurons

p = 0.1; % reservoir sparsity
k = 1;   % number of outputs

% input type
input_type = param.input_type;

%% Initialize matrix weights
% input weights
input = Win.*rand(N,N_th).*(rand(N,N_th) < Winp);

% static weights
rng('shuffle')
if Pexc == 0
    static = G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);
    
    % set the row average to be zero, explicitly, to induce chaotic spiking
    for i = 1:1:N
        QS = find(abs(static(i,:))>0);
        static(i,QS) = static(i,QS) - sum(static(i,QS))/length(QS);
    end
else
    % apply Dale's law
    disp('Dales law is applied')
    
    % number of excitatory and inhibitory neurons
    Nexc = round(Pexc * N);
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
end

% feedback weights
feedback = Q*(2*rand(N,k)-1);

% output weights
output = zeros(N, k);

% store all weights
weights.static = static;
weights.input = input;
weights.feedback = feedback;
weights.output = output;

%% List of trials for training and validation
train_trials = param.train_trials;
test_trials = param.test_trials;

%% Load input data
if param.makespikes || strcmp('ConvTrace', input_type) || strcmp('PSTH', input_type)
    % load the KernelStruct
    filename = ['.' f 'Input' f 'KernelStruct.mat'];
    
    if ~exist(filename, 'file')
        error('KernelStruct.mat is not in the input folder')
    end
    
    KernelStruct = load(filename);
    KernelStruct = KernelStruct.KernelStruct;
    
    % load the whiskmat
    filename = ['.' f 'Input' f 'whiskmat.mat'];
    
    if ~exist(filename, 'file')
        error('whiskmat.mat is not in the input folder')
    end
    
    whiskmat = load(filename);
    whiskmat = whiskmat.filtered_whiskmat;
end
%% Train and test the network
% run for N amount of epochs
for epoch = 1:N_total
    disp(['Epoch nr. ', num2str(epoch)])
    
    % ----------------------------TRAINING---------------------------------
    disp(['Training network, number of trials = ', num2str(N_train)])
    
    % apply FORCE learning
    FORCE = param.FORCE;
    weight_change = zeros(N_train, 1);
    
    for trial = 1:N_train
        
        % NETWORK INPUT
        trialId = train_trials(trial).trial;
        if param.makespikes || strcmp('ConvTrace', input_type) || strcmp('PSTH', input_type)
            
            % get the trial session and create the spikingstruct
            session = train_trials(trial).session;
            SpikeTrainStruct = make_trial_spikes(session, trialId,...
                whiskmat, KernelStruct);
        else
            % get the struct name and load it
            trial_mat = train_trials(trial).spike_struct;
            file = load( ['./Spiking structures/', trial_mat]);
            SpikeTrainStruct = file.SpikeTrainStruct;
        end
        
        % get the pole location and the input struct and target function
        pole = train_trials(trial).ytrain;
        [thalamus_input, target] =...
            reservoir_input(SpikeTrainStruct, input, N, N_th, pole, rate, input_type);
        
        % SIMULATE NETWORK
        % save the old output weights
        old_output = output;
        
        if strcmp('spikes', input_type)
            [~, output, ~, ~, ~] =...
                LIF_spiking_network(param, weights, thalamus_input, target, FORCE);
        else
            [~, output, ~, ~, ~] =...
                LIF_spiking_network_no_filt(param, weights, thalamus_input, target, FORCE);
        end
        
        % calculate the weight difference and update the output weights
        d_output = old_output - output;
        weight_change(trial, 1) = sum(abs(d_output));
        weights.output = output;
    end
    
    % --------------------------VALIDATION---------------------------------
    disp(['Testing network, number of trials = ', num2str(N_test)])
    
    % validation trials
    FORCE = false;
    for trial = 1:N_test
        
        % NETWORK INPUT
        trialId = test_trials(trial).trial;
        if param.makespikes || strcmp('ConvTrace', input_type) || strcmp('PSTH', input_type)
            
            % get the trial session and create the spikingstruct
            session = test_trials(trial).session;
            SpikeTrainStruct = make_trial_spikes(session, trialId,...
                whiskmat, KernelStruct);
        else
            % get the struct name and load it
            trial_mat = test_trials(trial).spike_struct;
            load( ['./Spiking structures/', trial_mat]);
        end
        
        % save the validation trials and first touches
        test_output.trials{trial} = trialId;
        test_output.first_touches(trial,1) = test_trials(trial).first_touch;
        
        % get the pole location and the input struct and target function
        pole = test_trials(trial).ytrain;
        [thalamus_input, target] =...
            reservoir_input(SpikeTrainStruct, input, N, N_th, pole, rate, input_type);
        
        % SIMULATE NETWORK
        if strcmp('spikes', input_type)
            [err, output_weights, Zx, Z_out, tspikes] =...
                LIF_spiking_network(param, weights, thalamus_input, target, FORCE);
        else
            [err, output_weights, Zx, Z_out, tspikes] =...
                LIF_spiking_network_no_filt(param, weights, thalamus_input, target, FORCE);
        end
        
        %TODO
        % input_save{trial}.neuron_input = thalamus_input;
        
        
        % get the spikinging statistics
        trial_length = length(Z_out);
        [A_t, ISI, Cv] = spike_stats(tspikes, trial_length , N);
        
        weights.output = output_weights;
        test_output.error{trial} = err;
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
    disp(['Test accuracy = ', num2str(acc)])
    training_output(epoch).weights = weights;
end

%% Save the output of the network
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
% saves the training output with all scaled parameter values as filename
savename = [savefolder f filename '.mat'];
save(savename, 'training_output', 'scale_param')

%TODO
% savename2 = [savefolder f 'input_save.mat'];
% save(savename2, 'input_save', '-v7.3');




