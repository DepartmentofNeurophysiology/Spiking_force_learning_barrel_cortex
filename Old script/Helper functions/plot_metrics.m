function  plot_metrics( perform, metrics ,alpha)

if metrics == 'val_acc'
    % Plot accuracy
    for a = 1:length(alpha)
        plot(perform(a).val_acc,'x-')
        hold on
    end
    hold off
    xlabel('Epoch')
    ylabel('Accuracy')
    if length(alpha) == 1
        legend(['Alpha = ' num2str(alpha(1))])
    elseif length(alpha) == 2
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))])
    elseif length(alpha) == 3
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))], ['Alpha = ' num2str(alpha(3))])
    end

elseif metrics == 'errors'
    % Plot errors for each trial
    for a = 1:length(alpha)
        plot(perform(a).errors,'x-')
        hold on
    end
    hold off
    xlabel('trial')
    ylabel('MSE')
    if length(alpha) == 1
        legend(['Alpha = ' num2str(alpha(1))])
    elseif length(alpha) == 2
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))])
    elseif length(alpha) == 3
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))], ['Alpha = ' num2str(alpha(3))])
    end
    
elseif metrics == 'sum_error'
    % Plot sum of errors for each epoch
    for a = 1:length(alpha)
        plot(perform(a).sum_error,'x-')
        hold on
    end
    hold off
    xlabel('Epoch')
    ylabel('Sum MSE')
    if length(alpha) == 1
        legend(['Alpha = ' num2str(alpha(1))])
    elseif length(alpha) == 2
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))])
    elseif length(alpha) == 3
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))], ['Alpha = ' num2str(alpha(3))])
    end
    
elseif metrics == 'weight_changes'
    % Plot absolute sum of weight changes for each epoch
    for a = 1:length(alpha)
        plot(perform(a).weight_changes,'x-')
        hold on
    end
    hold off
    xlabel('Epoch')
    ylabel('Sum MSE')
    if length(alpha) == 1
        legend(['Alpha = ' num2str(alpha(1))])
    elseif length(alpha) == 2
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))])
    elseif length(alpha) == 3
        legend(['Alpha = ' num2str(alpha(1))], ['Alpha = ' num2str(alpha(2))], ['Alpha = ' num2str(alpha(3))])
    end
end

end

