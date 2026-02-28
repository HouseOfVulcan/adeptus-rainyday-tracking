function p = pathFromRoot(varargin)
%pathFromRoot Build an absolute path under the configured Trackbench root.
%
% Usage:
%   p = trackbench.util.pathFromRoot("config", "default.json")

root = trackbench.util.rootDir();
p = fullfile(root, varargin{:});
end
