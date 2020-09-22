function [tspike, avg_fire_rate, Cv] = ni_nf_spiking_network(OMEGA, N)
% This function simulates a reservoir of spiking neurons with no input,
% output or feedback weights. 

% As input parameters it receives the static weights en the amount of neurons in the
% network. 
% As output returns the spike times, the average firing rates and the coefficient of variation. 

%% Network parameters
Ibias = -40; % current bias
dt = 0.05;% integration step time
tref = 2; % refractory time constant in milliseconds
tm = 10; % membrane time constant
vreset = -65; % voltage reset
vthresh = -40; % voltage threshold
td = 50; % synaptic decay
tr = 2; % synaptic rise

%% Define times and steps with the integration step time
T = 3500;
nt = T/dt;

%% Storage values

IPSC = zeros(N,1); %post synaptic current storage variable 
h = zeros(N,1); %Storage variable for filtered firing rates
r = zeros(N,1); %second storage variable for filtered rates 
hr = zeros(N,1); %Third variable for filtered rates 
JD = 0*IPSC; %storage variable required for each spike time  
v = vreset + rand(N,1)*(30-vreset); %Initialize neuronal voltage with random distribtuions
tlast = zeros(N,1); %This vector is used to set the refractory times
tspike = zeros(1,2);

%% MAIN NETWORK LOOP

% loop through the amount of network integration steps
for i = 1:1:nt
    
    % calculat the neuronal current
    I = IPSC + Ibias;
    
    % calculate the change in neuron voltage
    dv = (dt*i>tlast + tref).*(-v + I)/tm; 
    v = v + dt*(dv); 
    
    % remember the neurons that have spiked
    index = find(v>=vthresh); 
    
    % save voltage and current traces for every ms
    if mod(i,1/dt) == 0
        ms = i/(1/dt);
        vtrace(:,ms) = vtrace(:,ms) + v;
        Itrace(:,ms) = Itrace(:,ms) + I;
    end
   
    % if neurons have spiked
    if ~isempty(index)
        
        % calculate the sum synaptic input of the spikes
        JD = sum(OMEGA(:,index),2); 
        
        % save the spike times
        spikes = [index, 0*index+dt*i];
        tspike = [tspike, spikes];
    end
    
    % calculate the postsynaptic potential with double exponential filters
    IPSC = IPSC*exp(-dt/tr) + h*dt; 
    h = h*exp(-dt/td) + JD*(~isempty(index)/(tr*td));
    
    % calculate the output activity with double exponential filters
    %r = r*exp(-dt/tr) + hr*dt;
    %hr = hr*exp(-dt/td) + (v>=vthresh)/(tr*td);
    
    % reset the neurons that spiked and update the refactory period
    v = v + (30 - v).*(v>=vthresh)/(tr*td); 
    v = v + (vreset - v).*(v>=vthresh);
    tlast = tlast + (dt*i - tlast).*(v>=vthresh); 
end

%% Firing rate and coefficient of variation

% calculate the average firing rates
neurons = tspike(: , 1);
avg_fire_rate = calc_avg_fire_rate(neurons, N, T);

% calculate the coefficient of variation
Cv = calc_cv(tspike, N).';
end

