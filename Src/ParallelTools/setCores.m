%setCores will determine how many processors to turn on or off
%depending on the user specifications.
%
%  setCores(NumCores)
%
%  Success = setCores(NumCores)
%
%  [Success, ActiveCores] = setCores(NumCores)
%
%  INPUT
%    NumCores [N or 'max']: N or max number of processors to use.
%
%  OUTPUT
%    Success: 1 for successful or 0 for not successful
%    ActiveCores: Number of cores currently running

function [Success, ActiveCores] = setCores(NumCores)
MaxCores = getCores('max');
ActiveCores = getCores;
Success = 1;

if isempty(NumCores)
    return
end

if ischar(NumCores) 
    if strcmpi(NumCores, 'max')
        NumCores = MaxCores;
    else
        warning('%s: Unknown NumCores input "%s" - no changes made.', mfilename, NumCores);
        return
    end
end
    
NumCores = round(NumCores);
if NumCores > MaxCores
    NumCores = MaxCores;
elseif NumCores < 1
    NumCores = 1;
end

if ~isdeployed
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false; 
    ps.Pool.IdleTimeout = Inf;
end

if ActiveCores ~= NumCores 
    delete(gcp('nocreate'));
    if NumCores > 1
        try 
            parpool(NumCores);
            ActiveCores = NumCores;
        catch
            Success = 0;
        end
    end
end