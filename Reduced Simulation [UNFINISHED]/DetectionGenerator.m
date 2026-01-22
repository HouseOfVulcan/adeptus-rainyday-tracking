classdef DetectionGenerator < handle
    % =========================================================================
    % SENSOR SIMULATOR (RADAR/LIDAR)
    % =========================================================================
    %
    % PURPOSE:
    %   Acts as the interface between the Perfect World (Truth) and the 
    %   Noisy World (Tracker). It degrades ground truth data to simulate 
    %   hardware imperfections.
    %
    % PHYSICS SIMULATED:
    %   1. GEOMETRY: Targets must be within Range/Azimuth/Elevation limits.
    %   2. STOCHASTIC LOSS (Pd): Even visible targets are sometimes missed
    %      (e.g., Radar scintillation or scan gaps).
    %   3. NOISE (R): Positions are jittered by Gaussian noise defined by R.
    %   4. CLUTTER (False Alarms): "Ghost" dots are injected into the volume
    %      following a Poisson spatial distribution.
    %
    % =========================================================================
    
    properties
        Config  % Configuration Structure
    end
    
    methods
        function obj = DetectionGenerator(varargin)
            % --- DEFAULT CONFIGURATION ---
            defaultConfig.Pd             = 0.95;    % 95% Chance to see target
            defaultConfig.R              = diag([50, 50, 50].^2); % 50m Error
            defaultConfig.FalseAlarmRate = 1e-11;   % Density of clutter
            defaultConfig.MaxRange       = 50000;   % 50km
            defaultConfig.FOV            = [180; 90]; % [Az; El]
            defaultConfig.VolumeLimits   = [-20000 20000; -20000 20000; 0 10000];
            
            % Apply Overrides
            if nargin > 0 && isstruct(varargin{1})
                userConfig = varargin{1};
                fnames = fieldnames(userConfig);
                for i = 1:length(fnames)
                    defaultConfig.(fnames{i}) = userConfig.(fnames{i});
                end
            end
            obj.Config = defaultConfig;
        end
        
        function detections = step(obj, truth_list, time)
            detections = {};
            
            % 1. PROCESS REAL TARGETS
            for i = 1:length(truth_list)
                trueState = truth_list(i).State; 
                truePos   = [trueState(1); trueState(3); trueState(5)];
                
                % Visibility Check
                if ~obj.is_in_fov(truePos), continue; end
                
                % Probability of Detection Check (Miss Logic)
                if rand() > obj.Config.Pd, continue; end
                
                % Generate Gaussian Noise
                try L = chol(obj.Config.R, 'lower'); 
                catch, L = sqrt(obj.Config.R); end
                noise = L * randn(3, 1);
                
                det.Measurement      = truePos + noise;
                det.MeasurementNoise = obj.Config.R; 
                det.Time             = time;
                det.Type             = 'Target';
                detections{end+1} = det; %#ok<AGROW>
            end
            
            % 2. GENERATE CLUTTER (Poisson Point Process)
            volLimits = obj.Config.VolumeLimits;
            vol = diff(volLimits(1,:)) * diff(volLimits(2,:)) * diff(volLimits(3,:));
            numClutter = poissrnd(vol * obj.Config.FalseAlarmRate);
            
            for k = 1:numClutter
                clutX = volLimits(1,1) + rand() * diff(volLimits(1,:));
                clutY = volLimits(2,1) + rand() * diff(volLimits(2,:));
                clutZ = volLimits(3,1) + rand() * diff(volLimits(3,:));
                clutPos = [clutX; clutY; clutZ];
                
                if obj.is_in_fov(clutPos)
                    det.Measurement      = clutPos;
                    det.MeasurementNoise = obj.Config.R;
                    det.Time             = time;
                    det.Type             = 'Clutter';
                    detections{end+1} = det; %#ok<AGROW>
                end
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