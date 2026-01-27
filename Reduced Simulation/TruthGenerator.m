classdef TruthGenerator < handle
    % =========================================================================
    % TRUTH GENERATOR (FIXED: PERSISTENT TURBULENCE)
    % =========================================================================
    %
    % CHANGE LOG:
    %   - Added 'DriftState' to store accumulated wind effects.
    %   - Decoupled turbulence from Waypoints. Now, even if waypoints are 
    %     0.1s apart, the wind drift continues to push the aircraft off-course
    %     continuously over the full 60s.
    % =========================================================================
    
    properties
        Targets
        Config       % Weather Config (TurbulenceSigma)
        DriftStates  % Struct to store cumulative Pos/Vel drift per target
        LastTime     % To calculate dt for integration
    end
    
    methods
        function obj = TruthGenerator(weatherConfig)
            if nargin < 1
                obj.Config.TurbulenceSigma = 0; 
            else
                obj.Config = weatherConfig;
            end
            
            obj.Targets = struct('ID', {}, 'Waypoints', {}, ...
                                 'StartDelay', {}, 'ArrivalTimes', {}, ...
                                 'SegmentVels', {});
                                 
            % Initialize Drift Memory
            obj.DriftStates = containers.Map('KeyType', 'double', 'ValueType', 'any');
            obj.LastTime = -1; 
        end
        
        function add_target(obj, id, waypoints, speed_mps, start_delay)
            if nargin < 5, start_delay = 0; end
            
            nPoints = size(waypoints, 1);
            arrivalTimes = zeros(nPoints, 1);
            velocities   = zeros(nPoints-1, 3);
            arrivalTimes(1) = start_delay;
            
            for i = 1:nPoints-1
                p_curr = waypoints(i, :);
                p_next = waypoints(i+1, :);
                dist = norm(p_next - p_curr);
                dt_seg = dist / max(speed_mps, 1e-6);
                arrivalTimes(i+1) = arrivalTimes(i) + dt_seg;
                velocities(i, :) = (p_next - p_curr) / dt_seg;
            end
            
            newTgt.ID = id; newTgt.Waypoints = waypoints;
            newTgt.StartDelay = start_delay; newTgt.ArrivalTimes = arrivalTimes;
            newTgt.SegmentVels = velocities;
            obj.Targets(end+1) = newTgt;
            
            % Initialize Drift for this target (Pos=0, Vel=0)
            ds.Pos = [0, 0, 0];
            ds.Vel = [0, 0, 0];
            obj.DriftStates(id) = ds;
        end
        
        function current_states = get_states(obj, time)
            current_states = struct('ID', {}, 'State', {});
            
            % Calculate physics delta time
            if obj.LastTime < 0
                dt = 0; 
            else
                dt = time - obj.LastTime;
            end
            obj.LastTime = time;
            
            % Prevent huge jumps if time resets or lags
            if dt > 1.0, dt = 0; end 
            
            for i = 1:length(obj.Targets)
                tgt = obj.Targets(i);
                
                % Check if active
                if time >= tgt.ArrivalTimes(1) && time <= tgt.ArrivalTimes(end)
                    
                    % 1. CALCULATE IDEAL "RAIL" STATE
                    idx = find(tgt.ArrivalTimes <= time, 1, 'last');
                    if idx == length(tgt.ArrivalTimes)
                        ideal_pos = tgt.Waypoints(end, :); ideal_vel = [0, 0, 0];
                    else
                        t_start = tgt.ArrivalTimes(idx);
                        dt_seg  = time - t_start;
                        ideal_vel = tgt.SegmentVels(idx, :);
                        ideal_pos = tgt.Waypoints(idx, :) + (ideal_vel * dt_seg);
                    end
                    
                    % 2. UPDATE ACCUMULATED DRIFT (The Turbulence Fix)
                    ds = obj.DriftStates(tgt.ID);
                    
                    if dt > 0 && obj.Config.TurbulenceSigma > 0
                        % Random Acceleration (Gust)
                        accel = randn(1,3) * obj.Config.TurbulenceSigma;
                        
                        % Integrate: Accel -> Vel -> Pos
                        ds.Vel = ds.Vel + (accel * dt);
                        ds.Pos = ds.Pos + (ds.Vel * dt);
                        
                        % Simple Damping (Air Resistance)
                        % Prevents drift from accelerating to infinity over long runs
                        ds.Vel = ds.Vel * 0.99; 
                        
                        % Save State
                        obj.DriftStates(tgt.ID) = ds;
                    end
                    
                    % 3. COMBINE IDEAL + DRIFT
                    final_pos = ideal_pos + ds.Pos;
                    final_vel = ideal_vel + ds.Vel;
                    
                    s.ID = tgt.ID;
                    s.State = [final_pos(1); final_vel(1); final_pos(2); final_vel(2); final_pos(3); final_vel(3)];
                    current_states(end+1) = s;
                end
            end
        end
    end
end