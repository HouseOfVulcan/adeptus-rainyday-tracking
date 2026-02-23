function [results, detections] = runScenario(config, configName, detections)

    % RUNEXPERIMENT Execute one simulation run
    %
    % INPUT:
    %   config - Fully validated, normalized configuration struct
    %
    % OUTPUT:
    %   results    - Standard results struct (see ResultsSchema)
    %   detections - Detection log struct used by the trackers
    %
    % ASSUMPTIONS:
    %   - Config is valid (loader already checked)
    %   - Units are normalized (all meters, seconds, radians)
    %   - Coordinate frame is ENU
    %   - No plotting (caller handles visualization)
    %
    % what should the config look like
    % what is ENU
    %
    % call functions/classes that do things, don't do things
    % loader -> config -> runScnenario.m ->


    arguments
        % The fully loaded/validated configuration struct (required)
        config (1,1) struct

        % Optional name/identifier for the run
        configName (1,1) string = ""

        % Optional precomputed detections (to avoid regeneration)
        detections = []
    end
    
    %% Parameter Extraction and Setup
    % Extract Scenario Parameters
    numTargets    = config.scenario.num_targets;
    sceneDuration = config.scenario.duration_s;

    % Extract Degradation Setting
    enableDegradation = config.degradation.enabled;

    % Extract Global Tracker Parameters
    volume    = config.tracker_global.volume;
    beta      = config.tracker_global.beta;

    % Extract active tracker parameters (already selected by loadConfig)
    params = config.active_params;
    pd = params.pd;

    %% Generate Detections if not provided
    if isempty(detections)
        scenario = trackbench.scenario.createScenario( ...
            "NumTargets", numTargets, ...
            "SceneDuration", sceneDuration);

        detections = trackbench.detections.createDetections(scenario, enableDegradation);
    end

    %% Extract visualization flags with safe defaults
    if isfield(config.output, 'show_visuals')
        showVis = config.output.show_visuals;
    else
        showVis = true;
        fprintf("Warning: 'show_visuals' not specified in config.output, defaulting to true.\n");
    end

    if isfield(config.output, 'animate_visuals')
        animVis = config.output.animate_visuals;
    else
        animVis = true;
        fprintf("Warning: 'animate_visuals' not specified in config.output, defaulting to true.\n");
    end

    %% Visualization is handled by the entrypoint.

    %% Quick stats on detection count per scan
    % Helps confirm degradation is actually happening:
    %   IDEAL -> generally stable count close to number of targets
    %   RAINY -> increased variability (clutter + dropouts)
    if config.output.print_diagnostics
        nPerScan = cellfun(@numel, detections.Detections);
        fprintf("Detections/scan stats: min=%g, mean=%.2f, max=%g\n", ...
            min(nPerScan), mean(nPerScan), max(nPerScan));
    end

    fprintf(" | gate=%.1f | farGNN=%.2e | farMHT=%.2e | farJPDA=%.2e | pd=%.2f | volume=%.2e | beta=%.2e\n", ...
         params.gate, params.far_gnn, params.far_mht, params.far_jpda, pd, volume, beta);

    %% Run all enabled trackers dynamically
    results = trackbench.results.ResultsSchema.create();
    results.run_id = configName;
    results.config = config; % save the full config for this run in the results struct for traceability

    % Define the combinations we want to check based on the config
    trackerCombos = {
        {'GNN', 'CV',   config.trackers_to_run.gnn_cv};
        {'GNN', 'IMM',  config.trackers_to_run.gnn_imm};
        {'TOMHT','CV',  config.trackers_to_run.tomht_cv};
        {'TOMHT','IMM', config.trackers_to_run.tomht_imm};
        {'JPDA', 'CV',  config.trackers_to_run.jpda_cv};
        {'JPDA', 'IMM', config.trackers_to_run.jpda_imm}
        };

    for c = 1:length(trackerCombos)
        tType = trackerCombos{c}{1};
        fModel = trackerCombos{c}{2};
        isEnabled = trackerCombos{c}{3};

        if isEnabled
            comboName = lower(sprintf('%s_%s', tType, fModel));
            fprintf('\n============ %s + %s ============\n', tType, fModel);

            % 1. Use the new factory to build the tracker
            tracker = trackbench.tracking.buildTracker(tType, fModel, params, config.tracker_global, config.filter_params, pd);

            % 2. Run the tracker
            [trackSummary, truthSummary, trackMetrics, truthMetrics, time] = trackbench.tracking.runTracker(detections, tracker, false, showVis, animVis);

            % 3. Save results
            results.tracker_results.(comboName).trackSummary = trackSummary;
            results.tracker_results.(comboName).truthSummary = truthSummary;
            results.tracker_results.(comboName).trackMetrics = trackMetrics;
            results.tracker_results.(comboName).truthMetrics = truthMetrics;
            results.tracker_results.(comboName).time = time;

            if config.output.print_diagnostics
                disp(trackSummary); disp(truthSummary);
                disp(trackMetrics); disp(truthMetrics);
            end
        end
    end

    fprintf("\n==============================\n");
    fprintf(" RUN END\n");
    fprintf("==============================\n\n");
end
