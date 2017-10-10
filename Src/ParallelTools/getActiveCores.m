function ActiveCores = getActiveCores
Pool = gcp('nocreate');
if isempty(Pool)
    ActiveCores = 1;
else
    ActiveCores = Pool.NumWorkers;
end
