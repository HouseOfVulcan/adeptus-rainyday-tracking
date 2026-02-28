function config = loadConfig(configName, varargin)
%LOADCONFIG Backward-compatible wrapper around trackbench.config.loadConfig.
%   Keep legacy call sites working while centralizing all loader logic in
%   +trackbench/+config/loadConfig.m.

config = trackbench.config.loadConfig(configName, varargin{:});
end
