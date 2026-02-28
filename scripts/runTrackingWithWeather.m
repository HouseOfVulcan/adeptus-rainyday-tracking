function trackingWithWeather(configName)
% trackingWithWeather  Backward-compatible wrapper for runExperiment flow.
%
% This preserves the legacy script entrypoint while delegating execution to
% the package implementation.

arguments
    configName (1,1) string = "default_reRunDetections"
end

setupTrackbench();
trackbench.batch.runExperiment(configName);
end
