%findConsecLoc will take a binary array and locate the furthest X number of
%consecutive 1's away from the start location ('left' or 'right').
%
%  TrimLen = findConsecLoc(Q,MinConsec,StartSide) where Q is a binary 1 0
%  matrix, MinConsec is the number of consecutive matches to consider, and
%  StartSide is either 'left' or 'right'.
%
%  EX:
%    Q = [0 1 1 1 0 0 1 1 1 1 0 0];
%    MinConsec = 3;
%    TrimLen = findConsecLoc(Q,MinConsec,'left')
%      TrimLen = 10
%    TrimLen = findConsecLoc(Q,MinConsec,'right')
%      TrimLen = 11

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