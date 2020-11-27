function [run] = run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q,...
    Winp, alpha, Pexc, FORCE, makespikes, savefolder)
%RUN Summary of this function goes here
%   Detailed explanation goes here
%% FORCE training a reservoir of spiking neurons
% Explanation of the file here

%% Add subfolders and import data

addpath(genpath('Spiking structures'))
addpath(genpath('Helper data'))
addpath(genpath('Helper functions'))
addpath(genpath('Thalamus functions'))
addpath(genpath('Network functions'))
addpath(genpath('Parameter schemes'))

%% Parameters
tau_d = 50;             % synaptic decay
tau_r = 2;              % synaptic rise
Ibias = -40;            % bias current
step = 20;              % learning step
dt = 0.05;              % integration time constant
rate = 2;               % rate

%% Get the combination struct of the scaling parameters
param_comb = all_comb(Win, G, Q, Winp, Pexc);

%% Prepare a set of train and test trials

% load list of trial names
load('trainable_trials')

% get the shuffled train and test trials in the ratio 1:1, prox:dist
[train_trials, test_trials] = trial_selector(trainable_trials.prox_touch,...
    trainable_trials.dist_no_touch, N_train, N_test);

%TODO
%{
disp('fixed test trials, N_test = 2')
test_trials(1).trial = 80;  
test_trials(1).session = 'an171923_20120607';
test_trials(1).spike_struct = '101.mat';
test_trials(1).ytrain = [-1 77995];
test_trials(1).first_touch = 1892;
test_trials(1).pole_times = [1388 3342];

test_trials(2).trial = 129;
test_trials(2).session = 'an171923_20120613';
test_trials(2).spike_struct = '330.mat';
test_trials(2).ytrain = [1 152505];
test_trials(2).first_touch = 0;
test_trials(2).pole_times = [1354 3378];
%}
    
%% Prepare input struct
parameters_in = struct('N', N, 'N_th' , N_th, 'N_train', N_train,...
    'N_test', N_test, 'N_total', N_total, 'tau_d', tau_d,...
    'alpha', alpha, 'FORCE', FORCE, 'makespikes', makespikes,...
    'Ibias', Ibias, 'step', step, 'dt', dt, 'rate', rate, 'tau_r', tau_r,...
    'train_trials', train_trials, 'test_trials', test_trials);
%% Train and test the network

mkdir(savefolder)

parfor i = 1: size(param_comb, 1)
    
    % set the scaling parameters
    scale_parameters = struct();
    scale_parameters.Win = param_comb(i, 1);
    scale_parameters.G = param_comb(i, 2);
    scale_parameters.Q = param_comb(i, 3);
    scale_parameters.Winp = param_comb(i, 4);
    scale_parameters.Pexc = param_comb(i, 5);
    
    % run the network
    run(i).parameters_out = LIF_training(parameters_in,...
        scale_parameters, savefolder);
end

