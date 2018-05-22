%setPlotTickDecimal will extract the plot's current X and Y tick labels,
%and then set the number of decimal digits to be the same. This is useful
%since Matlab yields varying number of decimals points when plotting, which
%is not visually appealing.
%
%  setPlotTickDecimal(Ax, XDecimal, YDecimal)
%
%  setPlotTickDecimal(Ax, 'XDecimal', XDecimal, 'YDecimal', YDecimal)
%  
%  INPUT
%    Ax: axes handle to modify
%    Xdecimal [Int or 'max']: number of decimal digits only for the x axis.
%      Negative value will preserve current decimal digits. 'max' will find
%      the current maximum decimal count and apply it.
%    Ydecimal [Int or 'max']: number of decimal digits only for the y axis. 
%      Negative value will preserve current decimal digits. 'max' will find
%      the current maximum decimal count and apply it.

function varargout = setPlotTickDecimal(varargin)
[Axs, varargin] = getOpenGCA(varargin{:});
P = inputParser;
addParameter(P, 'Xdecimal', -1, @(x) isnumeric(x) || (ischar(x) && strcmpi(x, 'max')));
addParameter(P, 'Ydecimal', -1, @(x) isnumeric(x) || (ischar(x) && strcmpi(x, 'max')));

%Decide if user used the short (x, y) or long ('Xdecimal', x, ...) input.
if ~(~isempty(varargin) && ischar(varargin{1}) && ismember(lower(varargin{1}), {'getinput', 'getinputs', 'getvarargin'}))
    IsShortInput = 1;
    for j = 1:2:length(varargin)
        if ischar(varargin{j}) && ismember(varargin{j}, {'Xdecimal', 'Ydecimal', 'getinput', 'getinputs', 'getvarargin'})
            IsShortInput = 0;
            break;
        end
    end
    if IsShortInput %Make it the longer one
        NewVarargin = {'Xdecimal', -1, 'Ydecimal', -1};
        NewVarargin(2:2:length(varargin)*2) = varargin;
        varargin = NewVarargin;
    end
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInputT(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
Xdecimal = round(Ps.Xdecimal);
Ydecimal = round(Ps.Ydecimal);

%Use maximum decimal number if requested
if ischar(Xdecimal) && strcmpi(Xdecimal, 'max')
    Xdecimal = findMaxDecLen(get(Ax, 'XTickLabel'));
end
if ischar(Ydecimal) && strcmpi(Ydecimal, 'max')
    Ydecimal = findMaxDecLen(get(Ax, 'YTickLabel'));
end

%Do everything in inches to prevent rescale issues
PriorPositions = cell(length(Axs), 1);
for j = 1:length(Axs)
    PriorPositions{j} = get(Axs(j), 'Position');
end

%Set the X axis decimals
if Xdecimal >= 0
    for j = 1:length(Axs)
        %Decide if you want to use XTick value or XTickLabel value
        Xval = cellfun(@convStr2Num, get(Axs(j), 'XTickLabel'), 'unif', false);
        if any(cellfun(@(x) isempty(x) || isnan(x) || isinf(x), Xval))
            Xval = get(Axs(j), 'XTick');
        else
            Xval = cell2mat(Xval);
        end
        
        Xstr = cell(1, length(Xval));
        StrForm = ['%0.' num2str(Xdecimal) 'f'];
        for x = 1:length(Xstr)
            Xstr{x} = sprintf(StrForm, Xval(x));
        end
        set(Axs(j), 'XTickLabel', Xstr);
    end
end

%Set the Y axis decimals
if Ydecimal >= 0
    for j = 1:length(Axs)
        %Decide if you want to use YTick value or YTickLabel value
        Yval = cellfun(@convStr2Num, get(Axs(j), 'YTickLabel'), 'unif', false);
        if any(cellfun(@(x) isempty(x) || isnan(x) || isinf(x), Yval))
            Yval = get(Axs(j), 'YTick');
        else
            Yval = cell2mat(Yval);
        end
        
        Ystr = cell(1, length(Yval));
        StrForm = ['%0.' num2str(Ydecimal) 'f'];
        for y = 1:length(Ystr)
            Ystr{y} = sprintf(StrForm, Yval(y));
        end
        set(Axs(j), 'YTickLabel', Ystr);
    end
end

%Return the axes units, xlim, ylim to prior values
for j = 1:length(Axs)
    set(Axs(j), 'Position', PriorPositions{j});
end

%Finds the maximum number of deimals in a TickLabel
function MaxDec = findMaxDecLen(TickLabel)
MaxDec = 0;
for j = 1:length(TickLabel)
    DotLoc = find(TickLabel{j} == '.');
    if ~isempty(DotLoc)
        DecLen = length(TickLabel{j}) - DotLoc(end);
        if DecLen > MaxDec
            MaxDec = DecLen;
        end
    else
        continue
    end
end
