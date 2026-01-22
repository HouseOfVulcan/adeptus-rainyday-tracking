classdef TruthGenerator < handle
    % =========================================================================
    % TRUTH GENERATOR (KINEMATIC ENGINE)
    % =========================================================================
    %
    % PURPOSE:
    %   Generates deterministic Ground Truth states for targets moving between
    %   Waypoints at constant speed.
    %
    % IMPLEMENTATION:
    %   - Optimized "Look-up" architecture. Instead of integrating position 
    %     step-by-step, it pre-calculates the arrival time at every waypoint.
    %   - O(N) lookup during runtime to find current segment.
    %   - Outputs standard state vector: [x; vx; y; vy; z; vz].
    % =========================================================================
    
    properties
        Targets  % Struct Array of trajectory definitions
    end
    
    methods
        function obj = TruthGenerator()
            obj.Targets = struct('ID', {}, 'Waypoints', {}, ...
                                 'StartDelay', {}, 'ArrivalTimes', {}, ...
                                 'SegmentVels', {});
        end
        
        function add_target(obj, id, waypoints, speed_mps, start_delay)
            if nargin < 5, start_delay = 0; end
            
            % Pre-calculate segment timing
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
            
            newTgt.ID = id;
            newTgt.Waypoints = waypoints;
            newTgt.StartDelay = start_delay;
            newTgt.ArrivalTimes = arrivalTimes;
            newTgt.SegmentVels = velocities;
            obj.Targets(end+1) = newTgt;
        end
        
        function current_states = get_states(obj, time)
            current_states = struct('ID', {}, 'State', {});
            
            for i = 1:length(obj.Targets)
                tgt = obj.Targets(i);
                
                % Check if active
                if time >= tgt.ArrivalTimes(1) && time <= tgt.ArrivalTimes(end)
                    % Find current segment
                    idx = find(tgt.ArrivalTimes <= time, 1, 'last');
                    
                    if idx == length(tgt.ArrivalTimes)
                        pos = tgt.Waypoints(end, :);
                        vel = [0, 0, 0];
                    else
                        t_start = tgt.ArrivalTimes(idx);
                        dt_seg  = time - t_start;
                        vel = tgt.SegmentVels(idx, :);
                        pos = tgt.Waypoints(idx, :) + (vel * dt_seg);
                    end
                    
                    s.ID = tgt.ID;
                    s.State = [pos(1); vel(1); pos(2); vel(2); pos(3); vel(3)];
                    current_states(end+1) = s;
                end
            end
        end
    end
end