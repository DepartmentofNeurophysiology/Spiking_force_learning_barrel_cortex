function [ neuron_input, z_all ] = reservoir_input_trace( curve, angle, Ein, N, pole )

dat = pole;
whisk_input = [curve; angle];
%% Create the neuron input by multiplying with the input weights matrix

neuron_input = Ein*whisk_input;

%% Pulse parameters
reset = 1500; % Time inbetween the trials
pulse_length = 1000; % Length of pulse
amp = 2; % Amplitude of pulse
decay = 200; % decay of exponential pulse
constant = 500; % How long the pulse is kept constant

% if SpikeTrainStruct{1, 1}.first_touch == 0
%     first_touch = 1000;
% elseif SpikeTrainStruct{1, 1}.first_touch < 0
%     first_touch = 0;
% else
%     first_touch = SpikeTrainStruct{1, 1}.first_touch;
% end
first_touch = 500;

z_all = [];
start_early = 500; % How long to start before the end of the input
dt = 0.001;                                     % 1 msec                          
T_vec = 0:dt:reset*dt;                         % a vector with each time step	
%% Make the reservoir input and the target function

l_t = size(neuron_input,2) + reset;

z_t = zeros(l_t,1);
% z_t(size(spike_array.neuron_input,2)-start_early:1:size(spike_array.neuron_input,2)+constant) = dat{1}(1)*amp*ones(start_early+constant+1,1);
% z_t(size(spike_array.neuron_input,2)+constant:1:size(spike_array.neuron_input,2)+pulse_length+constant) =  dat{1}(1)*amp*exp(-(0:1:pulse_length)./decay);
z_t(first_touch:size(neuron_input,2)+constant) = -dat(1)*amp*ones(length(first_touch:size(neuron_input,2)+constant),1);
z_t(size(neuron_input,2)+constant:1:size(neuron_input,2)+pulse_length+constant) =  -dat(1)*amp*exp(-(0:1:pulse_length)./decay);

z_all = [z_all z_t'];

neuron_input = [neuron_input zeros(N, reset)];

end
