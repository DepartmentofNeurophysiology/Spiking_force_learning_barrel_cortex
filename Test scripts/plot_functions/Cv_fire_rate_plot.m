function Cv_fire_rate_plot(Cv , avg_fire_rate)
%CV_FIRE_RATE_PLOT Summary of this function goes here
%   Detailed explanation goes here

scatter(Cv, avg_fire_rate, 'filled')
xlabel('Cv')
ylabel('Firing rate Hz')

end

