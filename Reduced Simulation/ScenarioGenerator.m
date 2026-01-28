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
        function load_scenario1(truthGen)
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

        function load_scenario2(truthGen)
            % =============================================================
            % THE GNN STRESS TEST (Benchmark Scenario)
            % =============================================================
            % A comprehensive test suite for Adaptive GNN logic.
            %
            % TEST 1: THE JINK (Maneuver Response)
            %   - Target 1 flies straight, then executes a violent 90-degree 
            %     turn (High-G).
            %   - GOAL: Test if the 'LastNIS' logic triggers Q_Maneuver fast 
            %     enough to prevent track loss.
            %
            % TEST 2: THE WEAVE (Association Logic)
            %   - Targets 2 & 3 fly a tight sine-wave braid, crossing over 
            %     multiple times.
            %   - GOAL: Test if the 'Gate' and 'MatchPairs' logic can keep 
            %     ID 2 on T2 and ID 3 on T3 without swapping.
            %
            % TEST 3: THE SHADOW (Resolution/Merging)
            %   - Targets 4 & 5 fly parallel with only 100m separation 
            %     (inside the Merge Threshold), then split.
            %   - GOAL: Test if the reduced MergeThresholdSq allows them to 
            %     exist as separate tracks instead of getting merged.
            % =============================================================

            % --- TEST 1: THE JINK (Target 1) ---
            % 300 m/s. Flies Y-axis, Snaps to X-axis.
            wp_jink = [0, -10000, 5000;    % Start
                       0,      0, 5000;    % The Corner (Snap Turn)
                   10000,      0, 5000];   % Exit
            truthGen.add_target(1, wp_jink, 300, 0);

            % --- TEST 2: THE WEAVE (Targets 2 & 3) ---
            % High speed braid pattern
            wp_w1 = []; wp_w2 = [];
            amp = 800; freq = 0.001; speed_weave = 350;
            start_x = -5000; len = 15000;
            
            for x = start_x : 50 : (start_x + len)
                % Sinusoid Logic
                y_base = 2000;
                dy = amp * sin(freq * (x - start_x) * 2 * pi);
                
                wp_w1 = [wp_w1; x, y_base + dy, 5000];
                wp_w2 = [wp_w2; x, y_base - dy, 5000]; % 180 deg out of phase
            end
            truthGen.add_target(2, wp_w1, speed_weave, 0);
            truthGen.add_target(3, wp_w2, speed_weave, 0);
            
            % --- TEST 3: THE SHADOW (Targets 4 & 5) ---
            % Close formation flight (Parallel 80m separation)
            % This specifically tests your new "1.0 Sigma" merge fix.
            y_shadow = -4000; sep = 80;
            wp_s1 = [-10000, y_shadow, 4000; 
                       5000, y_shadow, 4000;   % Fly parallel until here
                       8000, y_shadow + 2000, 4000]; % Split Left
                       
            wp_s2 = [-10000, y_shadow + sep, 4000; 
                       5000, y_shadow + sep, 4000;   % Fly parallel until here
                       8000, y_shadow - 2000, 4000]; % Split Right
            
            truthGen.add_target(4, wp_s1, 280, 0);
            truthGen.add_target(5, wp_s2, 280, 0);
        end

        function load_scenario3(truthGen)
            % =============================================================
            % THE AIR SHOW (Showcase Scenario)
            % =============================================================
            % A visually complex scenario designed to look good in 3D plots
            % while testing 4 distinct tracking regimes simultaneously.
            %
            % 1. THE RIBBON (T1): A horizontal Figure-8.
            %    Tests: Continuous turning and self-intersection (crossing own path).
            %
            % 2. THE CORKSCREW (T2): A high-speed vertical spiral dive.
            %    Tests: 3D tracking and changing Z-velocity (accelerating dive).
            %
            % 3. THE BURST (T3-T5): 3 jets in tight formation that split apart.
            %    Tests: Resolution (merging) and multi-target initiation.
            %
            % 4. THE SCANNER (T6): A fast jet flying a sharp Sawtooth pattern.
            %    Tests: Rapid, instantaneous heading changes (High-G turns).
            % =============================================================

           % --- 1. THE RIBBON (T1: Figure-8) ---
            % RESCALED: Reduced size (a=3000) so it finishes in ~50s.
            wp_ribbon = [];
            a = 3000; % Reduced from 6000
            for t = 0 : 0.05 : 2*pi
                % Lemniscate Equations
                x = (a * cos(t)) / (1 + sin(t)^2);
                y = (a * sin(t) * cos(t)) / (1 + sin(t)^2);
                z = 5000 + (500 * sin(2*t)); 
                wp_ribbon = [wp_ribbon; x, y, z];
            end
            % Speed increased slightly to 350 m/s to ensure closure
            truthGen.add_target(1, wp_ribbon, 350, 0);

            % --- 2. THE CORKSCREW (Target 2: Spiral Dive) ---
            % Spirals down from 9000m to 1000m
            wp_screw = [];
            radius = 2000; center_x = -5000; center_y = 5000;
            for z = 9000 : -100 : 1000
                theta = (9000 - z) * 0.002; % Winding factor
                x = center_x + radius * cos(theta);
                y = center_y + radius * sin(theta);
                wp_screw = [wp_screw; x, y, z];
            end
            % Pull up at the end
            wp_screw = [wp_screw; center_x+5000, center_y+5000, 4000]; 
            truthGen.add_target(2, wp_screw, 350, 0);

            % --- 3. THE BURST (Targets 3, 4, 5: Formation Split) ---
            % Start: Tight V-Formation at (-10000, -5000)
            % End:   Split into 3 directions
            
            % Shared Start Path
            start_pt = [-10000, -5000, 3000];
            break_pt = [  0, -5000, 3000]; % Break point
            
            % T3 (Center -> Straight)
            wp_t3 = [start_pt; break_pt; 10000, -5000, 3000];
            
            % T4 (Left Wing -> Split Left/Up)
            wp_t4 = [start_pt + [0, 100, 0];  % 100m offset
                     break_pt + [0, 100, 0];
                     10000, -2000, 6000];     % Break Left/Up
            
            % T5 (Right Wing -> Split Right/Down)
            wp_t5 = [start_pt + [0, -100, 0]; % -100m offset
                     break_pt + [0, -100, 0];
                     10000, -8000, 1000];     % Break Right/Down
                     
            truthGen.add_target(3, wp_t3, 400, 0);
            truthGen.add_target(4, wp_t4, 400, 0);
            truthGen.add_target(5, wp_t5, 400, 0);

            % --- 4. THE SCANNER (T6: Combat Weave) ---
            % FIXED: Wavelength increased to 6000m (Physics realistic)
            wp_saw = [];
            y_base = 7000; z_base = 6000;
            amp = 2000;    
            
            for x = -9000 : 100 : 9000
                % Frequency = 1 / 6000m
                y = y_base + amp * sin(x * (1/6000) * 2 * pi);
                wp_saw = [wp_saw; x, y, z_base];
            end
            truthGen.add_target(6, wp_saw, 500, 0);

            % --- 5. THE MAVERICK (T7: The Thread) ---
            % A parabolic dive that swoops through the center of the scene.
            % Starts High (-12k, 9k), Dives to (0, 5k), Exits High (+12k, 9k).
            wp_mav = [];
            for x = -12000 : 100 : 12000
                % Normalized progress u (-1 to 1)
                u = x / 12000; 
                
                % Y: Sine wave to weave between zones (Amplitude 5000m)
                y = 5000 * sin(u * pi * 2); 
                
                % Z: Parabolic Dive (Starts 9000, bottoms at 4000, ends 9000)
                % At x=0, Z=4000 (Passes UNDER the Ribbon center)
                z = 4000 + 5000 * u^2; 
                
                wp_mav = [wp_mav; x, y, z];
            end
            truthGen.add_target(7, wp_mav, 500, 0);
        end
    end
end