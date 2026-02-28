function filter = custom_kalman(config)
% custom_kalman  User extension point for custom Kalman filter setup.

arguments
    config (1,1) struct
end

filter = struct();
filter.name = "custom_kalman";
filter.config = config;

error("custom_kalman:notImplemented", "Implement custom Kalman logic in extensions/trackers/custom_kalman.m");
end
