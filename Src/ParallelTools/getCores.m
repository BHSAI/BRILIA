%getCores will return the number of processor cores available or currently
%being used. By default, returns the currently active cores.
%
%  Cores = getCores;
%
%  Cores = getCores('max');
%
%  INPUT
%    'max': returns the maximum physical (not logical) cores available to
%      use on this system. Default is 'active', which returns the currently
%      active cores.
%
%  OUTPUT
%    Cores: the number of active (or maximum number of) cores
%
function Cores = getCores(varargin)
if contains('max', varargin', 'IgnoreCase', true)
    Cores = feature('numCores');
else
    Pool = gcp('nocreate');
    if isempty(Pool)
        Cores = 1;
    else
        Cores = Pool.NumWorkers;
    end
end