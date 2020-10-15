function [ dat ] = make_dat( s )
for i=1:length(s.trialIds)
    [dat(i).thetaVec, dat(i).kappaVec ,dat(i).timeVec, dat(i).pole,dat(i).touch_times,...
        dat(i).pole_vec,dat(i).pole_times,dat(i).lick_times ]...
     = angle_curve_position( s, i); % this creates structure for all trials with relevant info.
    if dat(i).pole_vec(1)==1 || dat(i).pole_vec(2)==1 % Label correct trials
        dat(i).correct = 1 ;        
    else
        dat(i).correct = 0 ;
    end
end

end

