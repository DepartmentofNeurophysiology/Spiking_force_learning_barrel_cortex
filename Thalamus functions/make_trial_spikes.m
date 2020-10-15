function [SpikingStruct] = make_trial_spikes(session, trialId)
%MAKE_TRIAL_SPIKES Summary of this function goes here
%   Detailed explanation goes here
f = filesep;

% load the KernelStruct
filename = ['.' f 'Input' f 'KernelStruct.mat'];

if ~exist(filename)
    error('KernelStruct.mat is not in the input folder')
end

KernelStruct = load(filename);
KernelStruct = KernelStruct.KernelStruct;

% load the whiskmat
filename = ['.' f 'Input' f 'whiskmat.mat'];

if ~exist(filename)
    error('whiskmat.mat is not in the input folder')
end

whiskmat = load(filename);
whiskmat = whiskmat.filtered_whiskmat;

% select sessions from the whiskingmat
session_index = find(strcmp({whiskmat.session}, session));
session_mat = whiskmat(session_index);

% select the trial from the sessions
trial_index = find([session_mat.trialId] == trialId);
trial_mat = session_mat(trial_index);

% create the SpikingStruct for specific trial
[SpikingStruct, ~] = dynamic_spike_maker(KernelStruct, trial_mat);
end

