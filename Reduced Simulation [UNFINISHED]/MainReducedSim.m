% MAINREDUCEDSIM.m
% =========================================================================
% TOP-LEVEL SIMULATION: "THE GRAND CHAOS"
% =========================================================================
% SCENARIO DESCRIPTION:
%   Combines "The Double Helix" and "The 5-Point Singularity".
%
%   1. THE HELIX (Targets 1 & 2): Two fighters perform a rolling scissors 
%      maneuver along the X-axis, intertwining constantly.
%
%   2. THE SINGULARITY (Targets 3-6): Four jets converge from North, South,
%      East, and West to meet at the exact center (0,0,5000).
%
%   3. THE CAMERA (Target 7): A drone spirals vertically around the
%      collision point to record the event.
%
%   CRITICAL MOMENT (T=28s): The Helix flies THROUGH the Singularity 
%   collision point just as the 4 jets merge.
% =========================================================================

clc; clear; close all;

%% 1. CONFIGURATION
sim.dt          = 0.25;   % High precision for the collision
sim.duration    = 60;     % Enough time for the full cross

% High-Fidelity Sensor settings
sensor.Pd             = 0.99;
sensor.FalseAlarmRate = 1e-12;   
sensor.R              = diag([100, 100, 100].^2); 

% Plot limits to see the whole 20km arena
plotConfig.Limits = [-10000 10000; -10000 10000; 0 10000];

%% 2. INITIALIZATION (THE GRAND CHAOS)
truthGen = TruthGenerator();

% -------------------------------------------------------------------------
% PART A: THE DOUBLE HELIX (Targets 1 & 2)
% -------------------------------------------------------------------------
% Path: Flies along the X-Axis from West to East
helix_len = 20000;  % 20km long
start_x   = -10000;
turns     = 5;
radius    = 1500;
center_z  = 5000;
speed_helix = 400;

wp_lead = []; wp_chase = [];

for x = start_x : 100 : (start_x + helix_len)
    prog = (x - start_x) / helix_len;
    theta = prog * (turns * 2 * pi);
    
    % Lead (Clockwise)
    y1 = radius * cos(theta);
    z1 = center_z + radius * sin(theta);
    
    % Chase (Counter-Clockwise offset)
    y2 = radius * cos(theta + pi); 
    z2 = center_z + radius * sin(theta + pi);
    
    wp_lead = [wp_lead; x, y1, z1];    %#ok<AGROW>
    wp_chase = [wp_chase; x, y2, z2];  %#ok<AGROW>
end

truthGen.add_target(1, wp_lead, speed_helix, 0);   % "Red Leader"
truthGen.add_target(2, wp_chase, speed_helix, 0);  % "Blue Two"

% -------------------------------------------------------------------------
% PART B: THE SINGULARITY (Targets 3, 4, 5, 6)
% -------------------------------------------------------------------------
% Four jets converging on (0,0,5000) from cardinal directions.
% They start 10km out and fly at 350 m/s.
% Impact time approx 28.5 seconds.

cross_dist = 10000;
speed_sing = 350;

% T3: East -> West (Intersects Helix Head-On)
truthGen.add_target(3, [cross_dist, 0, center_z; -cross_dist, 0, center_z], speed_sing, 0);

% T4: North -> South
truthGen.add_target(4, [0, cross_dist, center_z; 0, -cross_dist, center_z], speed_sing, 0);

% T5: South -> North
truthGen.add_target(5, [0, -cross_dist, center_z; 0, cross_dist, center_z], speed_sing, 0);

% T6: Vertical Dive (Top -> Down) for extra chaos
truthGen.add_target(6, [0, 0, 10000; 0, 0, 0], speed_sing, 0);

% -------------------------------------------------------------------------
% PART C: THE CAMERA DRONE (Target 7)
% -------------------------------------------------------------------------
% Spirals around the center to watch the crash
wp_cam = [];
cam_rad = 2500;
for ang = 0:10:720
    rad = deg2rad(ang);
    x = cam_rad * cos(rad);
    y = cam_rad * sin(rad);
    z = 2000 + (ang * 8); % Spirals up from 2km to 8km
    wp_cam = [wp_cam; x, y, z]; %#ok<AGROW>
end
truthGen.add_target(7, wp_cam, 300, 0);

% -------------------------------------------------------------------------
% COMPONENT SETUP
% -------------------------------------------------------------------------
detGen = DetectionGenerator(sensor);
tracker = GNNTracker(sim.dt); 
plotter = SimulationPlotter(plotConfig);
evaluator = PerformanceEvaluator();

%% 3. EXECUTION LOOP
fprintf('Starting "Grand Chaos" Simulation...\n');
numSteps = ceil(sim.duration / sim.dt);

for k = 1:numSteps
    time = (k-1) * sim.dt;
    
    % A. Generate Truth & Detections
    truth = truthGen.get_states(time);       
    dets  = detGen.step(truth, time);        
    
    % B. Update Tracker
    tracker.update(dets, sim.dt); 
    
    % C. Record & Plot
    evaluator.record_step(time, truth, tracker.Tracks);
    plotter.update(time, truth, dets, tracker.Tracks);
    
    % D. Console Status
    if mod(k, 10) == 0 
        nConf = sum(strcmp({tracker.Tracks.Status}, 'Confirmed'));
        fprintf('Time: %.1f | Truths: %d | Tracks: %d (Conf: %d)\n', ...
            time, length(truth), length(tracker.Tracks), nConf);
    end
    
    pause(0.01); 
end

fprintf('Simulation Complete.\n');

%% 4. FINAL REPORT
evaluator.evaluate(sim.dt);