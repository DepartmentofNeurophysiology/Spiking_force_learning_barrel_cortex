function [SpikingStruct] = make_trial_spikes(session, trialId, whiskmat, KernelStruct)
%MAKE_TRIAL_SPIKES Summary of this function goes here
%   Detailed explanation goes here

% select sessions from the whiskingmat
session_index = find(strcmp({whiskmat.session}, session));
session_mat = whiskmat(session_index);

% select the trial from the sessions
trial_index = find([session_mat.trialId] == trialId);
trial_mat = session_mat(trial_index);

% create the SpikingStruct for specific trial
[SpikingStruct, ~] = dynamic_spike_maker(KernelStruct, trial_mat);
end

