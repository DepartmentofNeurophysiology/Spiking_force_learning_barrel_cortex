function [changes ] = weight_change( parameters_out )

BPhi = parameters_out(1).BPhi{1,1};
changes = zeros(length(parameters_out)-1,1);
for i = 2: length(parameters_out)
    % Calculate difference in BPhi
    dBPhi = parameters_out(i).BPhi{1,1} - BPhi;
    % Set new BPhi for next iteration
    BPhi = parameters_out(i).BPhi{1,1};
    % calculate absolute sum of dBPhi and save
    changes(i-1) = sum(abs(dBPhi));
end


end

