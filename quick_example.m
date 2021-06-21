%% Quick example spiking force learning barrel cortex
% Runs short simulation of the network
%% Savefolder for output files
f = filesep;
savename = 'quick_example';
savefolder = ['.' f 'Output' f savename];
%% Set the parameters to run the simulation
% input parameters
N = 2000;        % number of neurons
N_th = 200;      % number of thalamus neurons
N_train = 2;     % number of training trials
N_test = 2;      % number of validation trials
N_total = 1;     % number of epochs

% weight scaling parameters
Win = 0.5;       % scales the input weights
G = 10;          % scales the static weights
Q = 1;           % scales the feedback weights
Winp = 1;        % network sparsity
alpha = 0.05;    % learning rate

% logicals
FORCE = true;       % if TRUE; apply FORCE learning during trials
makespikes = true;  % if TRUE; make the trial spiking structures 

% percentage of excitatory neurons if Dale's law is applied
Pexc = 0;   % percentage of excitatory neurons; set to 0 to ignore Dale's law restrictions

% Input 
input_type = 'spikes';  % options: 'ConvTrace', 'PSTH', 'spikes'
%% Run the simulation
run = run_sim(N, N_th, N_train, N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, input_type, savefolder);

