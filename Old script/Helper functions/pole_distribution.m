function   pole_distribution(s)
figure()
for i = 1:length(s.ytrain)
    if s.ytrain{i}(1) == -1
        plot(i,s.ytrain{i}(2),'b.') % proximal
        hold on
    else
        plot(i,s.ytrain{i}(2),'r.') % distal
        hold on
    end
end
xlabel('Trial number')
ylabel('Pole location')
legend('Proximal', 'Distal', 'Location','east' )
hold off
end

