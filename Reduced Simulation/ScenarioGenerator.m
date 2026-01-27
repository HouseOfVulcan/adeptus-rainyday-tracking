classdef ScenarioGenerator
    % SCENARIO GENERATOR
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