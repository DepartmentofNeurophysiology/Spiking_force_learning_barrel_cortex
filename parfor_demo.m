%% Varying parameters with parfor loops
% runs simulations of the network in a parfor loop with different 
% parameter values, saves run_sim output in a cell structure
% and calculates the computation time
%   * Varying number of training trials
%   * Varying number of epochs 
%   * Varying static weights

%% Add paths
f = filesep;
main_dir = 'parfor_demo_scaled';
mkdir(['Output' f main_dir])

%% Save files
savename = 'parfor';
savefolder_parfor = ['.' f 'Output' f main_dir f savename];

%% Global network parameters
% set the parameters for all simulations
% input parameters
N = 2000;       % number of neurons
N_th = 200;     % number of thalamus neurons
N_train = 500;  % number of training trials
N_test = 100;   % number of validation trials
N_total = 2;    % number of epochs 
alpha = 0.05;   % learning rate

% weight parameters
Win = 0.5;       % the input weights
G = 10;           % the static weights
Q = 1;           % the feedback weights
Winp = 1;        % network sparsity

% logicals
FORCE = true;       % apply FORCE learning during trials
makespikes = true;  % make the trial spiking structures 

% Dale's law
Pexc = 0;   % percentage of excitatory neurons; set to 0 to ignore Dale's law restrictions

% Input 
input_type = 'spikes';  % options: 'ConvTrace', 'PSTH', 'spikes'

%% Run the simulation
%% Varying number of training trials
run = cell(length(N_train));
N_train = 500:250:1500;
tic;
parfor i = 1:length(N_train)
    run{i} = run_sim(N, N_th, N_train(i), N_test, N_total, Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, input_type, savefolder_parfor)
end
toc;

%% Varying number of epochs
run = cell(length(N_total));
N_total = 1:5;
tic;
parfor i = 1:length(N_total)
    run{i} = run_sim(N, N_th, N_train, N_test, N_total(i), Win, G, Q, Winp,...
    alpha, Pexc, FORCE, makespikes, input_type, savefolder_parfor)
end
toc;

%% Varying static weights
run = cell(length(G));
G = 1:15;
tic;
parfor i = 1:length(G)
    run{i} = run_sim(N, N_th, N_train, N_test, N_total, Win, G(i), Q, Winp,...
    alpha, Pexc, FORCE, makespikes, input_type, savefolder_parfor)
end
toc;
