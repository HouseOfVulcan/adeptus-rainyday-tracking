% SCENARIOGENERATOR.m
% =========================================================================
% SCENARIO GENERATOR (Target Choreography)
% =========================================================================
% PURPOSE:
%   The "Director" of the simulation. It defines the exact flight paths 
%   (waypoints) for every target in the scene. By centralizing this logic, 
%   we guarantee that the "Single Run" debug mode and the "Monte Carlo" 
%   batch mode use the exact same physics and scenarios.
%
% HOW IT CONNECTS:
%   - Called by: MainReducedSim.m AND MonteCarloRunner.m.
%   - Modifies:  TruthGenerator (it calls 'add_target' to populate the world).
%
% THE SCRIPT (Default Scenario):
%   1. The Helix Pair (T1 & T2): Two fighters spiraling around each other. 
%      Tests the tracker's ability to handle high-speed, oscillating crossovers.
%   2. The Singularity (T3-T6): Four jets converging on a single point 
%      (0,0,5000) simultaneously. Tests the "Merge Logic" and collision handling.
%   3. The Camera Drone (T7): A slow, loitering observer flying a large 
%      circle. Acts as a "Control Group" (easy to track) to verify basic function.
% =========================================================================
% SCENARIO PARAMETER & TUNING GUIDE
% =========================================================================
% Use these variables to alter the geometry and difficulty of the test.
%
% 1. HELIX PARAMETERS (Agility Test)
%   - turns: Number of full rotations the targets make. Higher = More frequent crossings.
%   - radius: The width of the spiral (meters). Smaller radius = tighter turns = harder to track.
%   - speed_helix: Velocity (m/s). Higher speed requires a higher 'Q_Maneuver' in the tracker.
%
% 2. SINGULARITY PARAMETERS (Collision Test)
%   - cross_dist: The starting distance from the center. Determines how long 
%     the tracker has to lock on before the collision occurs.
%   - speed_sing: How fast they merge. 350 m/s is roughly Mach 1.
%
% 3. DRONE PARAMETERS (Baseline)
%   - cam_rad: Radius of the loiter circle. 
%   - ang step (10): Controls the smoothness of the circle.
% =========================================================================
classdef ScenarioGenerator
    % Centralized library of target trajectories.
    % Ensures MainReducedSim and MonteCarloRunner always use identical physics.
    
    methods(Static)
        function load_default(truthGen)
            % 1. HELIX PAIR (Targets 1 & 2)
            % High-speed crossing pattern with vertical oscillation
            helix_len = 20000; start_x = -10000; turns = 5; radius = 1500; center_z = 5000;
            speed_helix = 400;
            wp_lead = []; wp_chase = [];
            
            for x = start_x : 100 : (start_x + helix_len)
                prog = (x - start_x) / helix_len;
                theta = prog * (turns * 2 * pi);
                y1 = radius * cos(theta); z1 = center_z + radius * sin(theta);
                y2 = radius * cos(theta + pi); z2 = center_z + radius * sin(theta + pi);
                wp_lead = [wp_lead; x, y1, z1]; wp_chase = [wp_chase; x, y2, z2]; 
            end
            truthGen.add_target(1, wp_lead, speed_helix, 0);   
            truthGen.add_target(2, wp_chase, speed_helix, 0); 
            
            % 2. THE SINGULARITY (Targets 3, 4, 5, 6)
            % Four jets converging on a single point (0,0,5000) simultaneously
            cross_dist = 10000; speed_sing = 350; center_z = 5000;
            
            % Horizontal convergence
            truthGen.add_target(3, [cross_dist, 0, center_z; -cross_dist, 0, center_z], speed_sing, 0);
            truthGen.add_target(4, [0, cross_dist, center_z; 0, -cross_dist, center_z], speed_sing, 0);
            truthGen.add_target(5, [0, -cross_dist, center_z; 0, cross_dist, center_z], speed_sing, 0);
            
            % Vertical dive
            truthGen.add_target(6, [0, 0, 10000; 0, 0, 0], speed_sing, 0);
            
            % 3. SURVEILLANCE DRONE (Target 7)
            % Circular loiter pattern observing the chaos
            wp_cam = []; cam_rad = 2500;
            for ang = 0:10:720
                rad = deg2rad(ang);
                x = cam_rad * cos(rad); y = cam_rad * sin(rad); z = 2000 + (ang * 8); 
                wp_cam = [wp_cam; x, y, z]; 
            end
            truthGen.add_target(7, wp_cam, 300, 0);
        end
    end
end