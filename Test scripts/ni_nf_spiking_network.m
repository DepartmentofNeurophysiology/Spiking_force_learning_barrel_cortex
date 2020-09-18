function [output] = ni_nf_spiking_network(OMEGA, N)
%SPIKING_NETWORK_ Summary of this function goes here
%   Detailed explanation goes here

%% Network parameters
Ibias = -40; % current bias
dt = 0.05;% integration step time
tref = 1; % refractory time constant in milliseconds
tm = 10; % membrane time constant
vreset = -65; % voltage reset
vthresh = -40; % voltage threshold
td = 50; % synaptic decay
tr = 2; % synaptic rise

%% Define times

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
vtrace = zeros(N,T);
Itrace = zeros(N,T);

%% MAIN NETWORK LOOP


for i = 1:1:nt
    
    I = IPSC + Ibias; % neuronal current
    dv = (dt*i>tlast + tref).*(-v + I)/tm; %Voltage equation
    v = v + dt*(dv); 
    index = find(v>=vthresh); % find the neurons that have spiked
    
    % save voltage and current traces for every ms
    if mod(i,1/dt) == 0
        ms = i/(1/dt);
        vtrace(:,ms) = vtrace(:,ms) + v;
        Itrace(:,ms) = Itrace(:,ms) + I;
    end
   
    if ~isempty(index)
        JD = sum(OMEGA(:,index),2); % sum synaptic input of spikes
        
        % save spike times
        spikes = [index,0*index+dt*i];
        tspike = [tspike; spikes ];
    end
    
    tlast = tlast + (dt*i - tlast).*(v>=vthresh); % update the refactory periods
    
    IPSC = IPSC*exp(-dt/tr) + h*dt; % first synaptic filter
    h = h*exp(-dt/td) + JD*(~isempty(index)/(tr*td)); % second synaptic filter with decay time
    
    %r = r*exp(-dt/tr) + hr*dt;
    %hr = hr*exp(-dt/td) + (v>=vthresh)/(tr*td);
    
    v = v + (30 - v).*(v>=vthresh)/(tr*td); 
    v = v + (vreset - v).*(v>=vthresh); % reset the neurons that spiked
end

%% Firing rate and coefficient of variation

% calculate the average firing rates
neurons = tspike(: , 1);
avg_fire_rate = calc_avg_fire_rate(neurons, N, T);

% calculate the coefficient of variation
Cv = calc_cv(tspike, N).';

%% create output struct

output.tspike = tspike;
output.vtrace = vtrace;
output.Itrace = Itrace;
output.avg_fire_rate = avg_fire_rate;
output.Cv = Cv;
end

