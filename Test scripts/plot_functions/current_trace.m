function current_trace(Itrace, T)
%CURRENT_TRACE Summary of this function goes here
%   Detailed explanation goes here
x = 1 : T; 
plot(x, Itrace)
xlim([0 T])
legend
xlabel('t in ms')
ylabel('current pulse in nA')

end

