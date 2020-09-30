%% FORCE training a reservoir of spiking neurons
% Explanation of the file here

%% Add subfolders and import data

addpath(genpath('Spiking structures'))
addpath(genpath('Helper data'))
addpath(genpath('Helper functions'))
addpath(genpath('Network functions'))

% create an output folder
mkdir Output

%% Parameters

% set number of thalamus and reservoir neurons
N_th = 200;
N = 200;

% set the training, validation trials (has to be even)
N_train = 50;
N_test = 10;

% set the number of epochs
N_total = 1;

% set the scaling parameters of the input, static and feedback weights
Win = 0.5;
G = 8;
Q = 1; 

% Win sparsity
Winp = 1;

% set the synaptic decay and learning rate
tau_d = 50;
alpha = 0.05;

% apply FORCE learning
FORCE = 1;

% static parameters
Ibias = -40;
step = 20;
dt = 0.05;
rate = 3;
tau_r = 2;

%% Get the combination struct of the scaling parameters
param_comb = all_comb(Win, G, Q, Winp);

%% Prepare a set of train and test trials

% load list of trial names
load('trainable_trials')

% get the shuffled train and test trials in the ratio 1:1, prox:dist
[train_trials, test_trials] = trial_selector(trainable_trials.prox_touch,...
    trainable_trials.dist_no_touch, N_train, N_test);

%% Prepare input struct
parameters_in = struct('N', N, 'N_th' , N_th, 'N_train', N_train,...
    'N_test', N_test, 'N_total', N_total, 'tau_d', tau_d,...
    'alpha', alpha, 'FORCE', FORCE, 'Ibias', Ibias, 'step', step,...
    'dt', dt, 'rate', rate, 'tau_r', tau_r, 'train_trials', train_trials,...
    'test_trials', test_trials);
%% Train and test the network
for i = 1: size(param_comb, 1)
    
    % set the scaling parameters
    scale_parameters = struct();
    scale_parameters.Win = param_comb(i, 1);
    scale_parameters.G = param_comb(i, 2);
    scale_parameters.Q = param_comb(i, 3);
    scale_parameters.Winp = param_comb(i, 4);
    
    % run the network
    parameters_out(i) = LIF_training_v1(parameters_in, scale_parameters);
end

%% Save the output data

save('Output/new_trial_test', 'parameters_out')