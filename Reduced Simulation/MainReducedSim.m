% MAINREDUCEDSIM.m
% =========================================================================
% TOP-LEVEL SIMULATION: "THE GRAND CHAOS"
% =========================================================================

clc; clear; close all;

%% 1. USER CONFIGURATION
% Simulation Settings
sim.dt             = 0.50;   
sim.duration       = 60;     

% --- EXECUTION MODE ---
% Set TRUE for Batch Statistics (50x runs, no plots)
% Set FALSE for Single Debug Run (1x run, detailed plots)
sim.enableMonteCarlo = false; 
sim.numRuns          = 10;     

% SCENARIO: 'Clear', 'Rainy', 'Fog', 'Storm'
weatherType = 'rainy'; 

sensorConfig = WeatherGenerator.getParams(weatherType);

if sim.enableMonteCarlo
    % =====================================================================
    % MODE A: MONTE CARLO BATCH RUN
    % =====================================================================
    runner = MonteCarloRunner(sim, sensorConfig, weatherType);
    runner.run();
    
else
    % =====================================================================
    % MODE B: DETAILED SINGLE RUN
    % =====================================================================
    
    %% 2. INITIALIZATION
    plotConfig.Limits = [-10000 10000; -10000 10000; 0 10000];
    sim.enablePlotting = true;  

    truthGen = TruthGenerator(sensorConfig); 
    detGen   = DetectionGenerator(sensorConfig);
    tracker  = GNNTracker(sim.dt, sensorConfig.TrackerTuning); 
    evaluator = PerformanceEvaluator();

    % --- SCENARIO GENERATION (New Class) ---
    ScenarioGenerator.load_default(truthGen);

    if sim.enablePlotting
        plotter = SimulationPlotter(plotConfig);
    end

    %% 3. EXECUTION LOOP
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Starting Simulation: %s\n', weatherType); 
    fprintf('Progress: '); 

    numSteps = ceil(sim.duration / sim.dt);
    startTime = tic;        

    for k = 1:numSteps
        time = (k-1) * sim.dt;
        
        % A. Truth & Detections
        truth = truthGen.get_states(time);       
        dets  = detGen.step(truth, time, sim.dt); 
        
        % B. Tracker Update
        tracker.update(dets, sim.dt); 
        
        % C. Record History
        evaluator.record_step(time, truth, tracker.Tracks, dets);
        
        % D. Live Plotting
        if sim.enablePlotting
            plotter.update(time, truth, dets, tracker.Tracks);
        end
        
        % --- RETRO PROGRESS BAR ---
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

    %% 4. FINAL REPORT & PLOT
    evaluator.evaluate(sim.dt); 

    if ~sim.enablePlotting
        fprintf('Generating Post-Simulation Plot...\n');
        plotter = SimulationPlotter(plotConfig);
        plotter.plot_final_history(evaluator.History);
    end
end