function [avg_fire_rate] = calc_avg_fire_rate(neurons, N, T)
%CALC_FIRE_RATE Summary of this function goes here
%   Detailed explanation goes here

% define the edges of the bins
edges = 1 : (N + 1);

% sort the list of spikes and calc the firing rate
sorted_neurons = sort(neurons);
fire_rate = histcounts(sorted_neurons, edges);
avg_fire_rate = fire_rate / (T*10^-3);
end

