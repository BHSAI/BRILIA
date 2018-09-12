classdef ProgressTracker <  handle
    properties (Access = public)
        Iter
        MaxIter
        Mod
        Str
        DataQueue
        AfterEach
        StatusHandle
    end
    
    methods
        function O = ProgressTracker(MaxIter, Mod, Str, StatusHandle)
            if nargin < 2 || isempty(Mod)
                Mod = round(MaxIter/20);
            end
            if nargin < 3 || isempty(Str)
                Str = 'Progress';
            end
            if nargin < 4
                StatusHandle = [];
            end
            
            O.MaxIter = MaxIter;
            O.Mod = Mod;
            O.Str = Str;
            O.StatusHandle = StatusHandle;

            O.Iter = 0;
            O.DataQueue = parallel.pool.DataQueue;
            O.AfterEach = afterEach(O.DataQueue, @(~) O.incr);
        end
        
        function O = set(O, Prop, Value)
            O.(Prop) = Value;
        end
        
        function O = incr(O)
            O.Iter = O.Iter + 1;
            if mod(O.Iter, O.Mod) == 0 || O.Iter == O.MaxIter
                showStatus(sprintf('%s %d/%d', O.Str, O.Iter, O.MaxIter), O.StatusHandle);
            end
        end
        
        function delete(O)
            delete(O.AfterEach);
            delete(O.DataQueue);
        end
    end
end      