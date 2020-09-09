function [ thetaVec, kappaVec, timeVec, pole,touch_times,pos_vec,pole_times,lick_times ]...
    = angle_curve_position( s, trial_num)
% This function will take as input the structure and trial number and
% returns the corresponding whisker angle, curvature, pole location, and
% whisker touches. Also returns pole_vec so only correct trials can be
% selected.
%% some basic stuff
trialIdx = find(s.trialIds == trial_num);
if (length(trialIdx) == 0)
    disp('Requested trialNum not found; aborting.');
    return;
end
trialStartTime = s.trialStartTimes(trialIdx);

tsa = s.timeSeriesArrayHash.value{1};

%% get whisker angle, curvature, time vec
tsa = s.timeSeriesArrayHash.value{1};
ti = find(ismember(tsa.trial,trial_num));
if (length(ti) == 0)
    disp(['Warning: no whisker data for trial ', num2str(trial_num)]);
    timeVec = 0;
    thetaVec = 0;      % Whisker angle
    kappaVec = 0;       % Curvature
    
else
    timeVec = tsa.time(ti)-trialStartTime; 
    thetaVec = tsa.valueMatrix(1,ti);       % Whisker angle
    kappaVec = tsa.valueMatrix(2,ti);       % Curvature
end

%% pole position

pos_vec=s.trialTypeMat(:,trialIdx);         

if pos_vec(1) == 1 || pos_vec(3) == 1 || pos_vec(5) == 1 % Pole is left
    pole = 1;
else
    pole = -1;  % Pole is right
end
    
%% Touches
es = s.eventSeriesArrayHash.value{2};
ti = find(ismember(es.eventTrials{1},trial_num));
if (length(ti) > 0)
    touch_times = es.eventTimes{1}(ti) - trialStartTime;
else
    touch_times=0;
end
%% Bar in reach time
es = s.eventSeriesArrayHash.value{1};
ti = find(ismember(es.eventTrials,trial_num));
pole_times = es.eventTimes(ti) - trialStartTime;
%% Licks (decision made)
for bi=3:4
    es = s.eventSeriesArrayHash.value{bi};
    ti = find(ismember(es.eventTrials,trial_num));
    if (length(ti) > 0)
        lick_times = es.eventTimes(ti) - trialStartTime;
    else
        lick_times = 0;
    end
end

end

