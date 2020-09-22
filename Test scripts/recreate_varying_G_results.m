%% Recreate no input, no feedback spiking network results with varying G
% This script runs a spiking network without input, output or feedback with a
% varying static weights scaling parameter G. Thereafter it plots creates
% various plots to describe the network attributes

%% add file paths
addpath('plot_functions')
addpath('help_functions')

% create an output folder if it not yet exists
%if ~isfile('output')
    %mkdir output
%end
addpath('output')

%% Neurons and static weights

% number of neurons in the network
N = 2000;

% set the scalar array and sparsity
G = 0:1:20;
p = 0.1; 

% Create a struct with all different scaled static weight matrices
OMEGA_mat = zeros(N, N, length(G));

% initialize a the random weights to be scaled
rng('shuffle')
pre_OMEGA = (randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p);

% scale and save all the static weight matrices
for i = 1 : length(G)
    
    % initialize the random static weights
    
    OMEGA = G(i) * pre_OMEGA;
    
    % what is happening here 
    for j = 1:1:N 
        QS = find(abs(OMEGA(j,:))>0);
        OMEGA(j,QS) = OMEGA(j,QS) - sum(OMEGA(j,QS))/length(QS);
    end
    
    % save the weights in a matrix
    OMEGA_mat(:, :, i) = OMEGA;
end


%% Run the no input, no feedback network with alle the weight matrices

% storage struct
output_struct = {};

% do a parameter sweep
parfor i = 1 : length(G)
    disp(i)
    
    % choos the static weights
    static_weights = OMEGA_mat(:, :, i);
    
    % simulate the network
    [~, avg_fire_rate, Cv] = ni_nf_spiking_network(static_weights, N);
    
    % allocate space
    output_struct{i} = zeros(N, 2);
    
    % save the output 
    output_struct{i}(:, 1) = avg_fire_rate;
    output_struct{i}(:, 2) = Cv;
end

% save the ouput matrix
filename = 'output/output_struct';
save(filename)

% threedimensional storage matrix
%output_mat = zeros(N, 2, length(G));

% create an output matrix
%for i = 1 : length(G)
    %output_mat(:, :, i) = output_struct{i};
%end
    
%% prepare the data for plotting

% create an empty storage arrays
Cv = []; 
avg_fire_rate = []; 
Gvar = []; 
avg_Cv = zeros(length(G), 1);

% loop through the data files
for i = 1 : length(G)
    
    output_avg_fire_rate = output_struct{i}(:, 1)';
    output_Cv = output_struct{i}(:, 2)'; 
    
    % extract the CV and firing rate
    Cv = [Cv output_Cv];
    avg_fire_rate = [avg_fire_rate output_avg_fire_rate];
    Gvar = [Gvar (G(i) * ones(1, N))];
    avg_Cv(i) = mean(output_Cv);
    
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

%% plot log(G) against the CV
figure (2)
plot(log(G), avg_Cv)
title('CV relation to G')
xlabel('log G')
ylabel('average CV')

%% Histogram of CV neurons

% set the histogram edges
max_Cv = max(Cv);
edges = linspace(0, max_Cv, 100);

figure(4)

for i = 1 : length(G)
    single_Cv = output_struct{i}(:, 2);
    
    histogram(single_Cv, edges, 'Displaystyle', 'stairs', 'LineWidth', 3)
    title('Histogram of CV neurons')
    ylabel('# neurons')
    xlabel('CV')
    legend
    hold on
end

