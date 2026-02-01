%LOAD_CONFIG Load and validate a JSON config file into a MATLAB struct.
%
% Usage:
%   cfg = load_config("default");              % loads config/default.json
%   cfg = load_config("scenarios/clear_weather"); % loads config/scenarios/clear_weather.json
%
% OUTPUT
%   cfg : struct with validated config values

function cfg = load_config(jsonPath)

    arguments
        jsonPath (1,1) string
    end

    % Resolve path relative to project root
    if ~endsWith(jsonPath, ".json")
        jsonPath = jsonPath + ".json";
    end
    
    fullPath = fullfile(pwd, "config", jsonPath);
    
    % Check if file exists; if not, try fallback to default
    if ~isfile(fullPath)
        fprintf("[WARNING] Config file not found: %s\n", fullPath);
        fprintf("[INFO] Attempting to load default config...\n");
        fullPath = fullfile(pwd, "config", "default.json");
        
        if ~isfile(fullPath)
            error("Config file not found: %s. No default available.", fullPath);
        end
    end
    
    % Load and parse JSON
    fprintf("[CONFIG] Loading: %s\n", fullPath);
    rawText = fileread(fullPath);
    cfg = jsondecode(rawText);

    % --- Validate required structure ---
    if ~isfield(cfg, "scenario")
        error("Config missing required top-level key: 'scenario'");
    end
    s = cfg.scenario;

    required = ["mode", "num_targets", "duration_s"];
    for k = required
        if ~isfield(s, k)
            error("Config missing required key: scenario.%s", k);
        end
    end

    % --- Normalize types (jsondecode can produce char/double) ---
    cfg.scenario.mode = string(cfg.scenario.mode);
    cfg.scenario.num_targets = double(cfg.scenario.num_targets);
    cfg.scenario.duration_s = double(cfg.scenario.duration_s);

    % --- Simple sanity checks ---
    if ~(cfg.scenario.mode == "3D" || cfg.scenario.mode == "2D")
        error("scenario.mode must be '3D' or '2D'. Got: %s", cfg.scenario.mode);
    end
    if cfg.scenario.num_targets < 1 || mod(cfg.scenario.num_targets, 1) ~= 0
        error("scenario.num_targets must be a positive integer. Got: %g", cfg.scenario.num_targets);
    end
    if cfg.scenario.duration_s <= 0
        error("scenario.duration_s must be > 0. Got: %g", cfg.scenario.duration_s);
    end
    
    fprintf("[CONFIG] Validation passed.\n");
end