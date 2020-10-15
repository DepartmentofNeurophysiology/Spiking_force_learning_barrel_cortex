function   plot_network_target( epoch, VAL_BATCH, which_run,run )

network = [];
target = [];

%for n = 1:epochs
for i = 1:VAL_BATCH
    network = [network run(which_run).parameters_out(epoch ).network{i}'];
    target = [target run(which_run).parameters_out(epoch ).target{i}];
end
%end
plot(network)
hold on
plot(target)
hold off
xlabel('Time /ms')
ylabel('Output')
legend('network output', 'target')

end

