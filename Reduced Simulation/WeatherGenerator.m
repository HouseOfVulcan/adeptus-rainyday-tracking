classdef WeatherGenerator
    % WEATHER GENERATION (Updated for "Survivable Storm")
    
    methods(Static)
        function config = getParams(scenarioType)
            % --- BASELINES ---
            config.Pd = 0.99; 
            config.R = diag([50, 50, 50].^2); 
            config.FalseAlarmRate = 1e-12;
            config.TurbulenceSigma = 0; 
            config.AttenuationFactor = 0;    
            
            % --- DEFAULT TUNING ---
            tuning.Gate           = 45.0;     
            tuning.Q_Cruise       = 5.0;      
            tuning.Q_Maneuver     = 4000.0;   
            tuning.ManeuverThresh = 4.0;      
            tuning.MaxSpeed       = 1200.0;
            tuning.ConfirmTime        = 1.5; 
            tuning.PotentialCoastTime = 1.0; 
            
            switch lower(scenarioType)
                case 'clear'
                    tuning.ConfirmTime = 1.5; 
                    
                case 'rainy'
                    config.Pd = 0.98; 
                    config.FalseAlarmRate = 2e-10; 
                    config.R = diag([70, 70, 70].^2); 
                    config.AttenuationFactor = 2e-5; 
                    config.TurbulenceSigma = 5.0;    
                    
                    tuning.Gate = 25.0;        
                    tuning.Q_Cruise = 10.0;    
                    tuning.Q_Maneuver = 3000.0; 
                    tuning.ManeuverThresh = 3.0;
                    tuning.ConfirmTime = 2.0; 
                    tuning.PotentialCoastTime = 1.5; 
                    
                case 'fog'
                    config.Pd = 0.95; 
                    config.R = diag([60, 60, 60].^2);
                    config.FalseAlarmRate = 5e-11;
                    config.AttenuationFactor = 1e-5; 
                    config.TurbulenceSigma = 0;      
                    tuning.Gate = 40.0;
                    tuning.ConfirmTime = 2.5; 
                    
                case 'storm'
                    % --- SENSOR & PHYSICS ---
                    config.Pd = 0.90; 
                    config.R = diag([150, 150, 150].^2); 
                    
                    % KEEP REDUCED FAR: This is critical.
                    config.FalseAlarmRate = 1e-10;        
                    
                    config.AttenuationFactor = 5e-5; 
                    config.TurbulenceSigma = 20.0;   
                    
                    tuning.Gate = 25.0;        
                    tuning.Q_Cruise = 10.0;    
                    tuning.Q_Maneuver = 3000.0; 
                    tuning.ManeuverThresh = 3.0;
                    tuning.ConfirmTime = 2.0; 
                    tuning.PotentialCoastTime = 1.5;
                    
                otherwise
                    warning('Unknown weather type. Using Clear defaults.');
            end
            
            config.TrackerTuning = tuning;
           
        end
    end
end