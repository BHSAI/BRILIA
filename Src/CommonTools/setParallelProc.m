%setParallelProc will determine how many processors to turn on or off
%depending on the user specifications.
%
%  setParallelProc(NumProc)
%
%  Success = setParallelProc(NumProc)
%
%  INPUT
%    NumProc [N or 'max']: N number of processors to use, or 'max' number
%      of processors to use.
%
%  OUTPUT
%    Success: 1 for successful or 0 for not successful

function varargout = setParallelProc(NumProc)
Success = 0;
if isempty(NumProc) %Don't do anything if no input
    Success = 0;
else %Determine what to do
    %Ensuring the user input for NumProc is valid and in numbers
    if ischar(NumProc) %If user input 'max', set NumProc to max.
        if strcmpi(NumProc,'max')
            NumProc = feature('numCores');
        else
            warning('setParallelProc: Unknown NumProc input. Not changing processors.');
            Success = 0;
        end
    else %If user input a number, ensure it's a reasonable value
        NumProc = round(NumProc); %Must be integer
        if NumProc < 1; %Can't be less than 1
            warning('setParallelProc: NumProc must be >= 1. Setting to 1.');
            NumProc = 1; 
        end 
        if NumProc > feature('numCores'); %Can't be more than Max
            warning('setParallelProc: NumProc is greater than max cores. Setting to max cores.');
            NumProc = feature('numCores'); 
        end 
    end
    
    %Decide what to do about processing
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false; %Ensure that parfor is not automatically run.
    PoolName = gcp('nocreate');
    if isempty(PoolName) %No current pool, so open if needed
        if NumProc > 1
            try
                parpool(NumProc);
                Success = 1;
            catch
                Success = 0;
            end
        end
    else %Has a pool. check if NumProc differs from NumWorkers
        if PoolName.NumWorkers ~= NumProc %It does differ, so reconfigure
            delete(gcp('nocreate'));
            if NumProc > 1 %If parallel processing desired, open it.
                try 
                    parpool(NumProc);
                    Success = 1;
                catch
                    Success = 0;
                end
            end
        else %No need to change as it's already set
            Success = 1; 
        end
    end
end

if nargout >= 1
    varargout{1} = Success;
end
