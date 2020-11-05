function [ reservoir_input, target_function ] = reservoir_input( SpikeTrainStruct, n, Ein, N, dat, rate )
% Function that takes the nth thalamic spike times and makes neuron spiking
% input for the reservoir and returns this. 

% SpikeTrainStruct: structure containing PSTH, SpikeTimes and SpikeCount of
% a number of trials.
% n: selects the trials in the struct for which you want to make the
% reservoir input and the target function.
% Ein: is the input weights matrix.
% dat: contains the pole location.
% rate: contains the rate of the intermediate poisson firing. 

% output: is a structure containing the spiking reservoir input and the
% target function.

rng('shuffle')

%% Create the spiking array of 200 thalamus neurons for the nth trial
spike_array.trial = zeros(200, length(SpikeTrainStruct{1,1}.PSTH{1,n}));

for i = 1:200
    spike_array.trial(i,SpikeTrainStruct{1, 1}.SpikeTimes{i, n}) = ones(length(SpikeTrainStruct{1, 1}.SpikeTimes{i, n}) ,1);
end
%% Create the neuron input by multiplying with the input weights matrix
spike_array.neuron_input = zeros(N,length(spike_array.trial));
for t=1:length(spike_array.trial)
    index_input = find(spike_array.trial(:,t) == 1);  % Find thalamic input neurons that have spiked
    if length(index_input) ~= 0
        spike_array.neuron_input(:,t) = sum(Ein(:,index_input),2);
    end
end
%% Pulse parameters
reset = 1500; % Time inbetween the trials
pulse_length = 1000; % Length of pulse
amp = 2; % Amplitude of pulse
decay = 200; % decay of exponential pulse
constant = 100; % How long the pulse is kept constant

% if SpikeTrainStruct{1, 1}.first_touch == 0
%     first_touch = 1000;
% elseif SpikeTrainStruct{1, 1}.first_touch < 0
%     first_touch = 0;
% else
%     first_touch = SpikeTrainStruct{1, 1}.first_touch;
% end
first_touch = 500;

neuron_input = [];
z_all = [];
start_early = 500; % How long to start before the end of the input
dt = 0.001;                                     % 1 msec                          
T_vec = 0:dt:reset*dt;                         % a vector with each time step	
%% Make the reservoir input and the target function

l_t = size(spike_array.neuron_input,2) + reset;

z_t = zeros(l_t,1);
% z_t(size(spike_array.neuron_input,2)-start_early:1:size(spike_array.neuron_input,2)+constant) = dat{1}(1)*amp*ones(start_early+constant+1,1);
% z_t(size(spike_array.neuron_input,2)+constant:1:size(spike_array.neuron_input,2)+pulse_length+constant) =  dat{1}(1)*amp*exp(-(0:1:pulse_length)./decay);

z_t(first_touch:size(spike_array.neuron_input,2)+constant) = -dat(1)*amp*ones(length(first_touch:size(spike_array.neuron_input,2)+constant),1);
z_t(size(spike_array.neuron_input,2)+constant:1:size(spike_array.neuron_input,2)+pulse_length+constant) =  -dat(1)*amp*exp(-(0:1:pulse_length)./decay);

z_all = [z_all z_t'];

% make the poisson input
for n = 1:200
    vt = rand(size(T_vec) - [0 1]);
    spikes = (rate*dt) > vt;
    thalamus_poisson(n).spike_times = find( spikes == 1);
end
[ poisson_input] = make_poisson_spikes_weighted(N, Ein, reset, thalamus_poisson);

neuron_input = [neuron_input spike_array.neuron_input poisson_input];

reservoir_input = neuron_input;
target_function = z_all;


end

