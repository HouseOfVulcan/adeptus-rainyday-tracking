% DETECTIONGENERATOR.M
% =========================================================================
% DETECTION GENERATOR (Sensor Model)
% =========================================================================
% PURPOSE:
%   Simulates a Radar/Lidar sensor by transforming the "Truth" state into 
%   stochastic detections. It accounts for sensor limitations including 
%   Field of View (FOV), measurement noise, detection probability (Pd), 
%   and Poisson-distributed clutter (false alarms).
%
% HOW IT CONNECTS:
%   - Input:  Truth List (from TruthGenerator.m) containing exact positions.
%   - Output: Detection List (Cell Array) fed directly into GNNTracker.m.
%
% KEY LOGIC:
%   - Time-Invariant Clutter: Uses the 'dt' parameter to scale the clutter 
%     rate, ensuring that the number of false alarms per second remains 
%     constant regardless of the simulation's frame rate.
%   - Stochastic PD: Implements a probability of detection that can be 
%     optionally attenuated based on target range.
%   - Field of View: Filters both real targets and clutter points through 
%     azimuth, elevation, and range constraints.
% =========================================================================
% DETECTION GENERATOR PARAMETER & TUNING GUIDE
% =========================================================================
% This section defines the physical characteristics of the simulated 
% Radar/Lidar sensor. These parameters determine the "reliability" of the 
% data entering the tracker.
%
% 1. DETECTION PHYSICS
%   - Pd (Probability of Detection): The baseline chance (0.0 to 1.0) that 
%      a target is seen in a frame. Lower values simulate "missed blips."
%   - AttenuationFactor: Simulates signal degradation over distance 
%      (e.g., due to rain). Reduces Pd as the target moves further away.
%
% 2. CLUTTER (FALSE ALARMS)
%   - FalseAlarmRate: The density of false detections generated per 
%      unit volume, per second. Controls the amount of random noise/ghosts.
%      (e.g., 1e-11 is clear, 1e-9 is heavy storm).
%   - VolumeLimits: The 3D bounding box [min, max] for X, Y, Z coordinates
%      where clutter points are allowed to spawn.
%
% 3. MEASUREMENT UNCERTAINTY
%   - R (Measurement Noise Covariance): A 3x3 matrix defining the precision
%      of the sensor. Larger diagonal values mean the sensor is "noisier" 
%      and position readings will jitter more around the true location.
%
% 4. GEOMETRIC CONSTRAINTS
%   - MaxRange: The absolute maximum distance (meters) the sensor can reach.
%   - FOV (Field of View): [Azimuth; Elevation] limits in degrees. Targets 
%      outside this cone are strictly invisible, regardless of Pd.
% =========================================================================
classdef DetectionGenerator < handle

    properties
        Config
    end
    
    methods
        function obj = DetectionGenerator(userConfig)
            % Default values if not specified in initiazlization. 
            defaultConfig.Pd = 0.95; defaultConfig.R = diag([50, 50, 50].^2);
            defaultConfig.FalseAlarmRate = 1e-11; defaultConfig.MaxRange = 50000;
            defaultConfig.FOV = [180; 90]; defaultConfig.VolumeLimits = [-20000 20000; -20000 20000; 0 10000];
            defaultConfig.AttenuationFactor = 0; 
            
            if nargin > 0
                fnames = fieldnames(userConfig);
                for i = 1:length(fnames), defaultConfig.(fnames{i}) = userConfig.(fnames{i}); end
            end
            obj.Config = defaultConfig;
        end
        
        % Main randomizer
        function detections = step(obj, truth_list, time, dt)
            if nargin < 4, dt = 1.0; end % Fallback default
            
            volLimits = obj.Config.VolumeLimits;
            vol = diff(volLimits(1,:)) * diff(volLimits(2,:)) * diff(volLimits(3,:));
            
            % Scale clutter by dt. 
            % FAR is measured "Per Second", not "Per Scan".
            lambda = vol * obj.Config.FalseAlarmRate * dt;
            numClutter = poissrnd(lambda);
            
            maxDets = length(truth_list) + numClutter;
            emptyDet = struct('Measurement', [0;0;0], 'MeasurementNoise', obj.Config.R, ...
                              'Time', time, 'Type', '');
            detections = repmat(emptyDet, maxDets, 1);
            
            count = 0;
            
            % 1. PROCESS REAL TARGETS
            for i = 1:length(truth_list)
                trueState = truth_list(i).State; 
                truePos   = [trueState(1); trueState(3); trueState(5)];
                range     = norm(truePos);
                
                if ~obj.is_in_fov(truePos), continue; end
                
                currentPd = obj.Config.Pd * exp(-obj.Config.AttenuationFactor * range);
                if rand() > currentPd, continue; end
                
                try L = chol(obj.Config.R, 'lower'); catch, L = sqrt(obj.Config.R); end
                noise = L * randn(3, 1);
                
                count = count + 1;
                detections(count).Measurement = truePos + noise;
                detections(count).MeasurementNoise = obj.Config.R;
                detections(count).Time = time;
                detections(count).Type = 'Target';
            end
            
            % 2. PROCESS CLUTTER
            if numClutter > 0
                clutX = volLimits(1,1) + rand(numClutter, 1) * diff(volLimits(1,:));
                clutY = volLimits(2,1) + rand(numClutter, 1) * diff(volLimits(2,:));
                clutZ = volLimits(3,1) + rand(numClutter, 1) * diff(volLimits(3,:));
                
                for k = 1:numClutter
                    cPos = [clutX(k); clutY(k); clutZ(k)];
                    if obj.is_in_fov(cPos)
                        count = count + 1;
                        detections(count).Measurement = cPos;
                        detections(count).MeasurementNoise = obj.Config.R;
                        detections(count).Time = time;
                        detections(count).Type = 'Clutter';
                    end
                end
            end
            
            detections = detections(1:count);
            
            % Legacy Compatibility
            if count > 0
                detections = num2cell(detections);
            else
                detections = {};
            end
        end
        
        function visible = is_in_fov(obj, pos)
            range = norm(pos);
            if range > obj.Config.MaxRange, visible = false; return; end
            az = atan2d(pos(2), pos(1)); 
            el = asind(pos(3) / max(range, 1e-6)); 
            if abs(az) > obj.Config.FOV(1), visible = false; return; end
            if abs(el) > obj.Config.FOV(2), visible = false; return; end
            visible = true;
        end
    end
end