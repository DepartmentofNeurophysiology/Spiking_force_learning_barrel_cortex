function [ training_output ] = LIF_training_trace(param, scale_param, savefolder)
%LIF_TRAINING Summary of this function goes here
%   Detailed explanation goes here

f = filesep;

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

% set Pexc
Pexc = scale_param.Pexc;

% number of network outputs 
k = 1;

%% Initialize matrix weights

% input weights
input = Win.*rand(N,N_th).*(rand(N,N_th) < Winp); 

rng('shuffle')

% static weights

if Pexc == 0
    static =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);

    % set the row average to be zero, explicitly 
    for i = 1:1:N 
        QS = find(abs(static(i,:))>0);
        static(i,QS) = static(i,QS) - sum(static(i,QS))/length(QS);
    end
else
    % Apply Dale's law
    disp('Dales law is applied')
    
    % number of exci and ini neurons 
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

%% list of trials for training and validation
train_trials = param.train_trials;
test_trials = param.test_trials;

%% Load the whiskingmat for the traces
% load the whiskmat
filename = ['.' f 'Input' f 'whiskmat.mat'];

if ~exist(filename)
    error('whiskmat.mat is not in the input folder')
end

whiskmat = load(filename);
whiskmat = whiskmat.filtered_whiskmat;

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

        % make or load the spikes
        if param.makespikes
            disp('no spikes to be made')
        end
        
        % get the trial name and curve and angle traces
        session = train_trials(trial).session;
        trialId = train_trials(trial).trial;
        
        % select sessions from the whiskingmat
        session_index = strcmp({whiskmat.session}, session);
        session_mat = whiskmat(session_index);

        % select the trial from the sessions
        trial_index = [session_mat.trialId] == trialId;
        trial_mat = session_mat(trial_index);
        
        % get the pole and whiskertrace
        pole = train_trials(trial).ytrain;
        
        p = 1;
        [curve, angle] = make_whisker_trace(trial_mat, p);
 

        % get the pole location and the input struct and target function
        [neuron_input, target] =...
            reservoir_input_trace(curve, angle, input, pole); 
        
        % SIMULATE NETWORK
        % save the old output weights
        old_output = output;
        
        [~, output, ~, ~, ~] =...
            LIF_spiking_network_trace(param, weights, neuron_input, target, FORCE);
        
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

        % make or load the spikes
        if param.makespikes
            disp('no spikes to be made')
        end
        
        % get the trial name and curve and angle traces
        session = test_trials(trial).session;
        trialId = test_trials(trial).trial;
        test_output.first_touches(trial,1) = test_trials(trial).first_touch;
        
        % select sessions from the whiskingmat
        session_index = strcmp({whiskmat.session}, session);
        session_mat = whiskmat(session_index);

        % select the trial from the sessions
        trial_index = [session_mat.trialId] == trialId;
        trial_mat = session_mat(trial_index);
        
        angle = trial_mat.thetaVec;
        curve = trial_mat.kappaVec;

        % get the pole location and the input struct and target function
        pole = test_trials(trial).ytrain;
        [neuron_input, target] =...
            reservoir_input_trace(curve, angle, input, pole);
        
        % SIMULATE NETWORK
        [ err, output_weights, Zx, Z_out, tspikes ] =...
            LIF_spiking_network_trace(param, weights, neuron_input, target, FORCE);
        
        % get the spikinging statistics
        trial_length = length(Z_out);
        [ A_t, ISI, Cv ] = spike_stats( tspikes, trial_length , N );
        
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

