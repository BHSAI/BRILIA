%findConsecLoc will take a binary array and locate the "furthest" position
%of a consecutively matched segment that is away from either the left or
%right side of a logical array. Used mainly in trimGeneEdge.
%
%  TrimLen = findConsecLoc(Q,MinConsec,StartSide) 
%
%  INPUTS
%    Q: A 1xN logical matrix
%    MinConsec: Minimum number of consecutive matches to find
%    StartSide ['left','right']: Start looking from this side
%  
%  OUTPUTS
%    TrimLen: Number of nts to trim to cover the consective matched region.
%
%  EXAMPLE
%    Q = [0 1 1 1 0  1 1 1 0 1 1 1 1 0 0];
%    MinConsec = 3;
%    TrimLen = findConsecLoc(Q,MinConsec,'left')
%    TrimLen = 
%           13
%    TrimLen = findConsecLoc(Q,MinConsec,'right')
%    TrimLen = 
%           14
%
%  See also trimGeneEdge

function TrimLen = findConsecLoc(Q,MinConsec,StartSide)
if MinConsec <= 0
    error('Error: MinConsec must be a positive integer greater than 0');
end

%The one exception to the rule.
if sum(Q) == length(Q) && length(Q) == MinConsec
    TrimLen = length(Q);
    return
end

TrimLen = 0; 
c = 0; %Concsecutive counter

if MinConsec > 1
    if strcmpi(StartSide,'right')
        for k = length(Q)-1:-1:1
            if (Q(k+1) + Q(k)) == 2
                c = c+1;
                if c >= MinConsec
                    TrimLen = length(Q) - k + 1;
                end
            else
                c = 1;
            end       
        end
    elseif strcmpi(StartSide,'left')
        for k = 1:length(Q)-1
            if (Q(k+1) + Q(k)) == 2
                c = c+1;
                if c >= MinConsec
                    TrimLen = k + 1;
                end
            else
                c = 1;
            end       
        end
    else
        disp('Warning: Nothing has been done. Set StartSide correctly');
    end
else
    if strcmpi(StartSide,'right')
        for k = length(Q):-1:1
            if Q(k) == 1
                TrimLen = length(Q) - k + 1;
            end
        end
    elseif strcmpi(StartSide,'left')
        for k = 1:length(Q)
            if Q(k) == 1
                TrimLen = k;
            end       
        end
    end
end
