classdef DetectionGenerator < handle
    % DETECTION GENERATOR (Fixed: Time-Normalized Clutter)
    
    properties
        Config
    end
    
    methods
        function obj = DetectionGenerator(userConfig)
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
        
        % UPDATED ARGUMENTS: Added 'dt'
        function detections = step(obj, truth_list, time, dt)
            if nargin < 4, dt = 1.0; end % Fallback default
            
            volLimits = obj.Config.VolumeLimits;
            vol = diff(volLimits(1,:)) * diff(volLimits(2,:)) * diff(volLimits(3,:));
            
            % FIX: Scale clutter by dt. 
            % This converts FAR from "Per Scan" to "Per Second".
            % result: Smaller dt = Fewer clutter points per frame (Same total per second)
            lambda = vol * obj.Config.FalseAlarmRate * dt;
            numClutter = poissrnd(lambda);
            
            maxDets = length(truth_list) + numClutter;
            
            % Optimization: Struct Array
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