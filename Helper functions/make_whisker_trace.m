function [ Xnz_c, Xnz_a ] = make_whisker_trace( tr1, dat, pole)
% Selects the whisking curvature and angle when poles in reach an returns
% these two arrays.
% If pole is set to 1 then whisker trace is made for when pole is in reach.
% If pole is set to 0 then whisker trace is made for the whole trial.

% Whisker curvature
noNAN = dat(tr1).kappaVec;
noNAN(isnan(noNAN)) = 0;
trial_1 = noNAN;

pole_1 = dat(tr1).pole_times;
time_1 = dat(tr1).timeVec;

X_1 = zeros(round(time_1(end)),1);

i=1;
for t=1:length(X_1)
    X_1(t) = trial_1(i);
    if mod(t,2) == 0
        i = i + 1;
    end
    if i == length(trial_1)
        break
    end
end
if pole == 1
    % Only input between poles
    for t=1:length(X_1)
        if t < pole_1(1) || t > pole_1(2)
            X_1(t) = -10000;
        end
    end
    % Remove zeros
    Xnz_c = [];
    for t=1:length(X_1)
        if X_1(t) ~= -10000
            Xnz_c = [Xnz_c, X_1(t)];
        end
    end
else
    Xnz_c = X_1;
end
% Whisker angle
noNAN = dat(tr1).thetaVec;
noNAN(isnan(noNAN)) = 0;
trial_1 = noNAN;

pole_1 = dat(tr1).pole_times;
time_1 = dat(tr1).timeVec;
%

X_1 = zeros(round(time_1(end)),1);

i=1;
for t=1:length(X_1)
    X_1(t) = trial_1(i);
    if mod(t,2) == 0
        i = i + 1;
    end
    if i == length(trial_1)
        break
    end
end

% Only input between poles
if pole == 1
    for t=1:length(X_1)
        if t < pole_1(1) || t > pole_1(2)
            X_1(t) = -10000;
        end
    end
    % Remove zeros
    Xnz_a = [];
    for t=1:length(X_1)
        if X_1(t) ~= -10000
            Xnz_a = [Xnz_a, X_1(t)];
        end
    end
else
    Xnz_a = X_1;
end

Xnz_a = pi.*Xnz_a./180; % from degree to radians
end

