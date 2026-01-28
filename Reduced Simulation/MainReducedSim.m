% MAINREDUCEDSIM.m
% =========================================================================
% PROJECT FILE ARCHITECTURE:
% =========================================================================
%   - MainReducedSim.m     : The master controller script used to configur parameters and toggle 
%                            between single-run debugging and Monte Carlo analysis.
%   - MonteCarloRunner.m   : A batch execution engine that runs multiple simulations silently 
%                            to aggregate statistical performance metrics.
%   - ScenarioGenerator.m  : A centralized library that defines target trajectories 
%                            to ensure consistent physics across modes.
%   - WeatherGeneration.m  : A configuration factory that defines sensor noise, clutter density, 
%                            and turbulence levels for different environmental scenarios.
%   - TruthGenerator.m     : The physics engine that calculates true target states over time, 
%                            including the effects of accumulated wind drift.
%   - DetectionGenerator.m : A sensor model that transforms truth data into noisy radar 
%                            detections while generating time-normalized clutter points.
%   - GNNTracker.m         : The primary tracking logic that uses Global Nearest Neighbor 
%                            association and Kalman filtering to maintain trajectories.
%   - CubatureKalmanFilter.m: A specialized math library providing high-accuracy non-linear 
%                             prediction and correction steps for state estimation.
%   - PerformanceEvaluator.m: The analysis module that calculates industry-standard metrics 
%                             (MOTA, OSPA, RMSE) to grade tracker performance.
%   - SimulationPlotter.m  : The visualization engine responsible for generating live 
%                            3D plots and final 4-view history summaries.
% =========================================================================
% TOP-LEVEL SIMULATION CONTROLLER
% =========================================================================
% PURPOSE:
%   This is the master entry point for the Multi-Target Tracking simulation.
%   It orchestrates the generation of Truth trajectories,
%   Sensor Detections (Clutter, Noise), and the GNN Tracker logic.
%
% MODES OF OPERATION:
%   1. SINGLE RUN (Debug & Showcase Mode):
%      - Runs one simulation with detailed visualization enabled.
%      - Generates live 3D plots and a specific per-target performance report.
%      - Best for debugging tracker logic or visualizing specific failures.
%
%   2. MONTE CARLO (Batch & Performance Mode):
%      - Runs N simulations back-to-back without plotting.
%      - Aggregates statistics (MOTA, OSPA, Ghost Count) to prove system robustness.
%      - Best for final validation and generating report data.
%
% KEY USER PARAMETERS:
%   - dt:           Time step between frames (lower = more CPU, higher precision).
%   - weatherType:      Sets the difficulty (Noise, Clutter, Turbulence)
%                       and activates the preset tracker tunings. 
%   - enableMonteCarlo: The master switch between Debug (False) and Batch (True).
%   - enablePlotting: For single runs only. Animate objects real-time if true
% =========================================================================

clc; clear; close all;

%% 1. USER CONFIGURATION
% -------------------------------------------------------------------------
% SIMULATION TIMING
% -------------------------------------------------------------------------
sim.dt       = 0.1;   % [s] Time Step. 0.10s is high fidelity, 0.50s is standard radar.
sim.duration = 60;     % [s] Total duration of the scenario.

% -------------------------------------------------------------------------
% EXECUTION MODE
% -------------------------------------------------------------------------
% Set TRUE  -> Runs 'numRuns' times silently. Prints statistical averages.
% Set FALSE -> Runs 1 time with Live Plotting. Prints detailed breakdown.
sim.enableMonteCarlo = false; 
sim.numRuns          = 10;     % Number of iterations (Only used if enableMonteCarlo = true)

% -------------------------------------------------------------------------
% SCENARIO DIFFICULTY
% -------------------------------------------------------------------------
% Options: 'Clear' (Easy), 'Fog' (Medium), 'Rainy' (Hard), 'Storm' (Expert)
%   - 'Rainy' is the recommended "Hero" demo (Balanced Clutter/Survivability).
%   - 'Storm' will stress-test the tracker with high turbulence and false alarms.
weatherType = 'rainy'; 

% Load the specific sensor parameters (R, Pd, Clutter Rate) for this weather
sensorConfig = WeatherGenerator.getParams(weatherType);


%% 2. SIMULATION DISPATCHER
if sim.enableMonteCarlo
    % =====================================================================
    % MODE A: MONTE CARLO BATCH RUN
    % =====================================================================
    % Instantiates the Batch Runner to handle loop execution and stat aggregation.
    runner = MonteCarloRunner(sim, sensorConfig, weatherType);
    runner.run();
    
else
    % =====================================================================
    % MODE B: DETAILED SINGLE RUN
    % =====================================================================
    
    % --- INITIALIZATION ---
    plotConfig.Limits = [-10000 15000; -10000 15000; 0 10000];
    sim.enablePlotting = false;  % Force animated plotting ON for single run

    % Create System Objects
    truthGen = TruthGenerator(sensorConfig); 
    detGen   = DetectionGenerator(sensorConfig);
    tracker  = GNNTracker(sim.dt, sensorConfig.TrackerTuning); 
    evaluator = PerformanceEvaluator();

    % Load Standard Trajectories
    ScenarioGenerator.load_scenario3(truthGen);

    % Initialize Plotter if enabled
    if sim.enablePlotting
        plotter = SimulationPlotter(plotConfig);
    end

    % --- EXECUTION LOOP ---
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Starting Simulation: %s\n', weatherType); 
    fprintf('Progress: '); 

    numSteps = ceil(sim.duration / sim.dt);
    startTime = tic;        

    for k = 1:numSteps
        time = (k-1) * sim.dt;
        
        % A. GENERATE DATA
        truth = truthGen.get_states(time);       
        dets  = detGen.step(truth, time, sim.dt); 
        
        % B. UPDATE TRACKER
        tracker.update(dets, sim.dt); 
        
        % C. LOG DATA
        evaluator.record_step(time, truth, tracker.Tracks, dets);
        
        % D. VISUALIZE (animate if on)
        if sim.enablePlotting
            plotter.update(time, truth, dets, tracker.Tracks);
        end
        
        % Progress Bar (Only show dots if plotting is off to avoid console spam)
        if ~sim.enablePlotting
            if mod(k, ceil(numSteps/50)) == 0
                fprintf('.'); 
            end
        end
    end

    totalTime = toc(startTime); 

    fprintf(' Done.\n'); 
    fprintf('\nSimulation Complete.\n');
    fprintf('Total Simulation Run Time: %.4f seconds\n', totalTime);
    fprintf('Real-Time Factor: %.2fx (Sim Time / Run Time)\n', sim.duration / totalTime);

    % --- FINAL REPORT ---
    evaluator.evaluate(sim.dt); 

    % If animation off, generate a static summary plot at the end
    if ~sim.enablePlotting
        fprintf('Generating Post-Simulation Plot...\n');
        plotter = SimulationPlotter(plotConfig);
        plotter.plot_final_history(evaluator.History);
    end
end