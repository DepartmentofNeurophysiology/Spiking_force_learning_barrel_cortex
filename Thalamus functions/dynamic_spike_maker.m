function [ SpikeTrainStruct, pole ] = dynamic_spike_maker(KernelStruct, dat)
% DYNAMIC_SPIKE_MAKER makes thalamic spikes trains for specific trials
% Input:
%   * KernelStruct = struct containing thalamic kernels
%   * dat = struct containing the wisking data
% Output:
%   * SpikeTrainStruct = struct containing the corresponding thalamic
% spike trains
%   * pole = pole (-1 or 1)

%% Make barrelstruct
Nbx = 1;                                                % # of barrels 'x-direction'
Nby = 1;                                                % # of barrels 'y-direction'
barrelstruct = cell(Nbx, Nby);
for nbx = 1:Nbx
    for nby = 1:Nby
        barrelstruct{nbx, nby}.xpos         = (nbx-2)*300;      % barrel position
        barrelstruct{nbx, nby}.ypos         = (nby-1)*300;
        barrelstruct{nbx, nby}.Nthalamic    = 200;              % # filter neurons for this barrel
        barrelstruct{nbx, nby}.mainbarrel   = 3;
    end
end
barrelstruct{1,1}.mainbarrel    = 1; % main
SpikeTrainStruct = cell(Nbx,Nby);
SpikeGenStruct = cell(Nbx,Nby);
nb = 0;
pole = 1; 
%% Get the whisker traces from the data
[WhiskerTrace.Recording{2,1} ,WhiskerTrace.Recording{1,1}]  = make_whisker_trace(dat, pole);
%WhiskerTrace.Recording{2,1} = dat.x_c{i};
%WhiskerTrace.Recording{1,1} = dat.x_a{i};
%% Make the thalamic spikes
%%
WhiskerTrace.binsize = 1;
Barrelstruct = barrelstruct;
for nbx = 1:Nbx
    for nby = 1:Nby
        nb = nb+1;
        SpikeGenStruct{nbx,nby}.refra             = 3; % refractory period (ms)
        SpikeGenStruct{nbx,nby}.Ntrial_pertrace   = 1; % # trials for each Deflection trace
        SpikeGenStruct{nbx,nby}.binsize           = 1; % binsize spike trains (ms)
        if Barrelstruct{nbx,nby}.mainbarrel == 1
            
            SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = 1; % Scaling of PSTH
        elseif Barrelstruct{nbx,nby}.mainbarrel == 2
            
            SpikeGenStruct{nbx,nby}.delay             = 2.5; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = .3;
        elseif Barrelstruct{nbx,nby}.mainbarrel == 3
            
            SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
            SpikeGenStruct{nbx,nby}.scaling           = 0;
        end
        
        plotyn = 1;
        SpikeTrainStruct{nbx,nby} = kernel_recording_to_spiketrain(WhiskerTrace, KernelStruct{nbx,nby}, SpikeGenStruct{nbx,nby}, [], plotyn);
    end
end

pole = dat.pole;
%pole = dat.ytrain(trials,1);
end

