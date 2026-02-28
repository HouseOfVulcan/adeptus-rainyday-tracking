function root = rootDir(varargin)
%rootDir Get/set Trackbench project root directory.
%
% Usage:
%   root = trackbench.util.rootDir()
%   root = trackbench.util.rootDir("/abs/path/to/trackbench")

key = "trackbench_root_dir";

if nargin > 0
    candidate = char(string(varargin{1}));
    if ~isfolder(candidate)
        error('trackbench:rootDir:notFound', 'Root directory does not exist: %s', candidate);
    end
    setappdata(0, key, candidate);
end

if isappdata(0, key)
    root = getappdata(0, key);
    return;
end

% Fallback auto-detect from this file location.
thisFile = mfilename('fullpath');
root = fileparts(thisFile);

for i = 1:8
    if isfile(fullfile(root, "config", "default.json"))
        setappdata(0, key, root);
        return;
    end
    parent = fileparts(root);
    if strcmp(parent, root)
        break;
    end
    root = parent;
end

% Final fallback preserves previous behavior for edge cases.
root = pwd;
setappdata(0, key, root);
end
