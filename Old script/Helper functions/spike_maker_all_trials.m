function  spike_maker_all_trials( dat, KernelStruct )
% Make spike for all trials
for i = 1:length(dat.trial_number)
    [ SpikeTrainStruct, pole ] = dynamic_spike_maker( i, KernelStruct, dat);
    SpikeTrainStruct{1,1}.ytrain = pole;
    SpikeTrainStruct{1,1}.struct_name = dat.struct_name(i);
    SpikeTrainStruct{1,1}.trial_number = dat.trial_number(i);
    save(['Spiking structures\',num2str(i)],'SpikeTrainStruct')
end

end

