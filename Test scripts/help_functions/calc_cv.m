function [Cv] = calc_cv(tspike, N)
%CALC_CV Summary of this function goes here
%   Detailed explanation goes here

% calculate the ISI
ISI = [];

% find the difference in spike time between the spikes for each neuron
for n = 1 : N
    index = tspike(:, 1) == n;
    spike_neuron(n).time = tspike(index, 2);
    spike_neuron(n).ISI = diff(spike_neuron(n).time);
    ISI = [ISI diff(spike_neuron(n).time)'];
end

% calculate the coefficient of variation for each neuron
Cv = zeros(N,1);

for n = 1:N
    Cv(n,1) = std(spike_neuron(n).ISI)/mean(spike_neuron(n).ISI);
end
end

