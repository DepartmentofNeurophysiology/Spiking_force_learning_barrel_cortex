%% Add all subfolders and import data.
addpath(genpath('Network loop functions'))
addpath(genpath('Helper functions'))
addpath(genpath('Thalamus functions'))
addpath(genpath('Spiking structures')) % THIS IS THE FOLDER WITH THE SPIKING STRUCTURES
addpath(genpath('Helper data'))
addpath(genpath('Input data'))
%% Set number of neurons
N_th = 200;                               % This is the number of thalamus neurons
N = 400;                                  % Number of neurons in the reservoir.
%% Parametersweep over these

G = [7, 8]; % scales static weights
Q =[1]; % scales learnt/output weights
par_comb = allcomb(G, Q);

%% Set other parameters

Win = [ 0.5 ];  % Scale the input weights.
Ein = Win.*rand(N,N_th); % Make input matrix

td = [50];% synaptic decay
alpha= [0.05]; % Learning rate


N_train = 2; % Number of trials to train on ( Has to be even)
N_test = 2; % Number of trials to test on ( Has to be even)
N_total = 2; % Number of epochs
%% Make one set of train and test ( validation) trials
load('trainable_trials') % List of trial names

prox = prox_touch(randperm(length(prox_touch))); % shuffle proximal trials
dist = dist_no_touch(randperm(length(dist_no_touch))); % shuffle distal trials
train_nr = N_train;
train_prox = prox(1:train_nr/2); % Make sure that 50 % of trials is distal and other is proximal
train_dist = dist(1:train_nr/2);
train_all = [train_prox train_dist];
val_prox = prox( train_nr/2 +1 : train_nr/2+ N_test/2);
val_dist = dist( train_nr/2 +1 : train_nr/2+ N_test/2);
val = [val_prox val_dist];

shuffled_val = val(randperm(length(val))) ; % Shuffle 
shuffled_train = train_all(randperm(length(train_all))) ; % Shuffle 

%% Define all input paramters in one input structure

parameters_in = struct('par_comb',par_comb,'N',N,'alpha',alpha,'BIAS',-40,'step',20,'dt',0.05,'rate',5, 'G', G, 'Q', Q, 'td', td, 'Ein', Ein,'shuffled_val',shuffled_val,'shuffled_train',shuffled_train);
parameters_in.N_train = N_train; % Number of trial to train on.
parameters_in.N_test = N_test; % Number of trials to test on for val. acc. Should be even
parameters_in.N_total = N_total; % Number times N_train and N_test is repeated

%% Train and test the network
foldername = 'G_8_4'; % To save data

%%rmdir(foldername)
mkdir(foldername)
clear parameters_out

% to run one loop at a time -> (i = 1: size(par_comb,1), 0)
parfor i = 1: size(par_comb,1) % Use parfor loop to train networks with different parameters
    
    run(i).parameters_out  = LIF_train_some( parameters_in, i, foldername); % All relevant parameters are stored in the 'run' struct
    
end

%% Save and sort all of the output data
%send_mail('Simulatie klaar!') % Send email to 

run_acc = zeros(length(run), length(run(1).parameters_out)); % To store the test accuracies in

%% This loop selects the info you want to save from the run struct. 
% As the run struct tends to get very large, MATLAB does not like saving
% this. Therefore this sections selects important data from the struct
% which is then saved. 

for i =1:length(run)
    for v  = 1:length(run(i).parameters_out)
        run_acc(i,v) = run(i).parameters_out(v).val_acc; % The test accuracy 
        dw{i,v} = run(i).parameters_out(v).weight_change; % Change in weights
        A_t{i,v} = run(i).parameters_out(v).A_t; % Activity of the network
        Cv{i,v} = run(i).parameters_out(v).Cv; % Coefficient of variation
        weights{i,v} = run(i).parameters_out(v).weights; % Weights of the network
    end
end

% Save this data to the foldername you have defined

save(foldername, 'run_acc', 'dw', 'A_t', 'Cv', 'parameters_in', 'weights')

%% The following section saves the network output, target, and reservoir spikes.

epoch = parameters_in.N_total;

for r =  1:size(par_comb,1)
    network = run(r).parameters_out(epoch).network;
    target = run(r).parameters_out(epoch).target;
    spikes = run(r).parameters_out(epoch).spikes;
    val_trials = run(r).parameters_out(epoch).val_trials;
    save_name = ['output/', num2str(r)];
    %save_name_target = ['target_run_', num2str(r)];
    save(save_name, 'network', 'target', 'spikes','val_trials')
end




%% Define trial you want to investigate
parameters_out = run(1).parameters_out;
epoch = 2;
trial = 2;
tspike = parameters_out(epoch).spikes{trial,1};A_t = parameters_out(epoch).A_t{trial,1};ISI = parameters_out(epoch).ISI{trial,1};Cv = parameters_out(epoch).Cv{trial,1};val_trial = parameters_out(epoch).val_trials{1,trial};  
Cv = Cv(~isnan(Cv));

%% Plot the spikes of the reservoir
figure(1)
subplot(2,1,1)
plot(tspike(:,2),tspike(:,1),'k.')
ylabel('Neuron')
xlabel('time/ ms')
ylim([0 200])
title('Network spikes')

subplot(2,1,2)
network = parameters_out(epoch).network{trial,1};
target = parameters_out(epoch).target{trial,1};
plot(network)
hold on
plot(target)
hold off
%%
