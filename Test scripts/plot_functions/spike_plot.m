function spike_plot(neurons, spike_times, T)
%Function that creates a spike plot from a sequence of neuron spikes and
%their spike times.
%   takes an array of a sequence of neuron spikes and their spike times and
%   the measurement time T as input

% create the plot
plot(spike_times, neurons,'k.')

% set the plot axis limits and labels
xlim([0 T])
title('Network Spikes')
xlabel('T in ms')
ylabel('neurons')

end

