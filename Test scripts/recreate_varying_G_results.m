%% Recreate no input, no feedback spiking network results with varying G
% This script runs a spiking network without input or feedback with a
% varying static weights scaling parameter (G) 

%% add file paths

addpath('plot_functions')
addpath('help_functions')
addpath('output')

%% Neurons and static weights

% number of neurons in the network
N = 2000;

% Create a struct with all different scaled static weight matrices (1-100)
% set the scalar array and sparsity
G = 1 : 100;
p = 0.1; 

% define the struct and weight matrix before scaling
OMEGA_struct = {};

% scale and save all the static weight matrices
for i = 1 : length(G)
    
    rng('shuffle')
    pre_OMEGA = (randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);
    
    OMEGA = G(i) * pre_OMEGA;
    
    for j = 1:1:N % What is happening here?
        QS = find(abs(OMEGA(j,:))>0);
        OMEGA(j,QS) = OMEGA(j,QS) - sum(OMEGA(j,QS))/length(QS);
    end
    
    OMEGA_struct{i} = OMEGA;
end


%% Run the no input, no feedback network with alle the weight matrices

output = {};

parfor i = 1 : length(OMEGA_struct)
    
    static_weights = OMEGA_struct{i};
    
    output{i} = ni_nf_spiking_network(static_weights, N);
    
    disp(i)
end

%% Save the data

for i = 1 : length(output)
    single_output = output{i};
    
     % save the output in the output file
    filename = strcat('output/', int2str(i), '.mat');
    save(filename, 'single_output');
end
    
%% plot the data

% create an empty storage arrays
Cv = []; %zeros(N*100, 1);
avg_fire_rate = []; % zeros(N*100, 1);
Gvar = [];  %zeros(N*100, 1);

% set the step
index = 1;
step = 10;

% loop through the data files
for i = 1:1:100
    
    file = strcat('output/', int2str(i), '.mat');
    output = load(file);
    
    % extract the CV and firing rate
    Cv = [Cv output.single_output.Cv];
    avg_fire_rate = [avg_fire_rate output.single_output.avg_fire_rate];
    Gvar = [Gvar (i * ones(1, 2000))];
    
    
end

%% plot the Cv against the firing rate over different values of G


figure (1)
scatter3(Cv, avg_fire_rate, Gvar, 15, Gvar, 'filled')
colormap(cool)
colorbar
title('No input, no feedback')
xlabel('CV')
xlim([0 6])
ylabel('firing rate in Hz')
ylim([0 40])

