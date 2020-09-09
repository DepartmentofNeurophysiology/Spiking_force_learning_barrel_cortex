function [ parameters_out ] = LIF_train_some( param, iter, foldername)
f = filesep;
% param: strcuture containing all the important parameters for the
% SNN.
% iter: tells spiking_network.m function which parameters to use. 
% foldername is to save data after one epoch, in case something goes wrong
% after this point.
%% Define paramters, some are fixed and some vary depending on the parameter search

% It is important that in this section the correct parameters are defined. 

% THIS DEPENDS ON THE PARAMETER SEARCH THAT YOU WANT TO PERFORM

% In this case these parameters vary per network
G = param.par_comb(iter,1); Q = param.par_comb(iter, 2); 

% And these are identical for each network
N = param.N;   
Ein = param.Ein;
rate = param.rate;
N_train = param.N_train; N_total = param.N_total; N_test = param.N_test;
k = 1; % number of outputs

%% Initialise weights
p = 0.1; %Set the network sparsity 
BPhi = zeros(N,k); % Output weights

rng('shuffle')
E = (2*rand(N,k)-1)*Q;  % Feedback weights
OMEGA =  G*(randn(N,N)).*(rand(N,N)<p)/(sqrt(N)*p); % Static weights
for i = 1:1:N
    QS = find(abs(OMEGA(i,:))>0);
    OMEGA(i,QS) = OMEGA(i,QS) - sum(OMEGA(i,QS))/length(QS);
end

%% Save weights to structure which will be input to the network
weights.output_weights = BPhi;
weights.static_weights = OMEGA;
weights.feedback_weights = E;
weights.input_weights = Ein;
%% Get all the file names of the spiking structures
filePattern = fullfile('..', f ,'Spiking structures', '*.mat'); 
theFiles = dir(filePattern);
names = load('labeled_spike_names.mat'); % These are the file names and their y and touch values.
names = names.names;
%names = names(:,randperm(length(names))); % Shuffle the stuff
%shuffled_names = names(randperm(length(names))) ; % Shuffle the file names
%% Make train set and validation set. This was already made in the script.

shuffled_val = param.shuffled_val;
shuffled_train = param.shuffled_train;

%% Run network on train set with RLS to update output weights & test network on validation set to find classification accuracy.

% TRAINING PART
for n = 1:N_total %Repeat everyting N_total times 
    disp(['EPOCH = ', num2str(n)])
    weight_change = zeros(N_train,1);
    for trial = 1:N_train
 
        trial_train = names{1, shuffled_train(trial)};
        load( ['.', f , 'Spiking structures' , f ,trial_train]) % Load the spiking structure from the spiking directory
        pole = SpikeTrainStruct{1,1}.ytrain;
        [ input_struct ] = reservoir_input( SpikeTrainStruct, 1, Ein, N, pole, rate ); % Get spiking input 
        FORCE = true; % RLS weight updates
        
        % MAIN NETWORK FUNCTION
        
        BPhi_old = BPhi;
        [~,BPhi,~,~,~] = spiking_network( param, weights, input_struct, FORCE, iter);
        dBPhi = BPhi_old - BPhi;
        weight_change(trial,1) = sum(abs(dBPhi));    % Sum of absolute weight change values.
        weights.output_weights = BPhi;               % Update the output weights
        
    end
    
    % Validate trials and save the accuracy
    
    FORCE = false; % No RLS weight updates, only look at the output of the network.
    
% VALIDATION PART

    clear stats
    clear val_dat
    save_val_trials = {};
    first_touches = zeros(N_test, 1);
    
    for trial = 1: N_test 

        trial_val = names{1,shuffled_val(trial)};
        save_val_trials{trial} = trial_val;
        load( ['.' , f , 'Spiking structures' , f ,trial_val]) % Load the spiking structure from the spiking directory
        first_touches(trial,1) = SpikeTrainStruct{1, 1}.first_touch;
        pole = SpikeTrainStruct{1,1}.ytrain;
        % Get spiking input and target funtion for a specfic trial
        [ input_struct ] = reservoir_input( SpikeTrainStruct, 1, Ein, N, pole, rate );
        
        % MAIN NETWORK FUNCTION
       
        [val_dat{trial,1}, val_dat{trial,2}, val_dat{trial,3}, val_dat{trial,4}, val_dat{trial,5}] = spiking_network( param, weights, input_struct, FORCE, iter);      
        
        [stats{trial,1}, stats{trial,2}, stats{trial,3}] = spike_stats( val_dat{trial,5}, length(val_dat{trial,4}) , N );
    end
    
    % Save data to structure and run again    

    parameters_out(n).error = val_dat(:,1);
    parameters_out(n).weights = weights;
    parameters_out(n).target = val_dat(:,3);
    parameters_out(n).network = val_dat(:,4);
    parameters_out(n).spikes = val_dat(:,5);
    parameters_out(n).weight_change = weight_change;
    parameters_out(n).A_t= stats(:,1);
    parameters_out(n).ISI= stats(:,2);
    parameters_out(n).Cv= stats(:,3);
    parameters_out(n).val_trials = save_val_trials;
    parameters_out(n).train_trials = shuffled_train;
    [ acc ] = val_acc( parameters_out(n), first_touches);
    
    parameters_out(n).val_acc = acc;
    
    % Save some stuff every epoch in case it stops somewhere
    
    clear save_struct
    save_struct = {'CV', stats(:,3), 'acc', acc};
    save_name = [foldername , f , 'run_' , num2str(iter) ,'epoch_' , num2str(n)];
    save(save_name, 'save_struct')
    
end
end

