%setPlotTickDecimal will extract the plot's current X and Y tick labels,
%and then set the number of decimal digits to be the same. This is useful
%since Matlab yields varying number of decimals points when plotting, which
%is not visually appealing.
%
%  setPlotTickDecimal(Gx, XDecimal, YDecimal)
%
%  setPlotTickDecimal(Gx, Param, Value, ...)
%  
%  INPUT
%    Gx: figure or axes handle
%    XDecimal [int or 'max']: number of decimal digits only for the x axis.
%      Negative value will preserve current decimal digits. 'max' will find
%      the current maximum decimal count and apply it.
%    YDecimal [int or 'max']: number of decimal digits only for the y axis. 
%      Negative value will preserve current decimal digits. 'max' will find
%      the current maximum decimal count and apply it.
%
%    Parameter    Values(*)   Description
%    -----------  ----------  ---------------------------------------------
%    XDecimal     * -1        Do not modify decimals in x axis
%                   int       Number of decimals to have in the x axis
%                   'max'     Use maximum decimals of current x axis values
%    YDecmial     * -1        Do not modify decimal in y axis
%                   int       Number of decimals to have in the y axis
%                   'max'     Use maximum decimals of current x axis values
%    ZDecmial     * -1        Do not modify decimal in z axis
%                   int       Number of decimals to have in the z axis
%                   'max'     Use maximum decimals of current x axis values
%
%  EXAMPLE
%    Gx = loadDemoPlot
%    setPlotTickDecimal(Gx, 4, 3)
%
function varargout = setPlotTickDecimal(varargin)
[Axs, varargin] = getOpenGCA(varargin{:});
P = inputParser;
addParameter(P, 'XDecimal', -1, @(x) isnumeric(x) || (ischar(x) && strcmpi(x, 'max')));
addParameter(P, 'YDecimal', -1, @(x) isnumeric(x) || (ischar(x) && strcmpi(x, 'max')));
addParameter(P, 'ZDecimal', -1, @(x) isnumeric(x) || (ischar(x) && strcmpi(x, 'max')));
if ~any(cellfun('isclass', varargin, 'char')) || any(strcmpi(varargin, 'max'))%Using short input style
    NewVarargin = {'XDecimal', -1, 'YDecimal', -1, 'ZDecimal', -1};
    NewVarargin(2:2:length(varargin)*2) = varargin;
    varargin = NewVarargin;
end
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return
end

PriorPositions = cell(length(Axs), 1);
for j = 1:length(Axs)
    PriorPositions{j} = get(Axs(j), 'Position');
end
D = {'X', 'Y', 'Z'};
for j = 1:length(Axs)
    for d = 1:length(D)
        if ischar(Ps.([D{d} 'Decimal'])) && strcmpi(Ps.([D{d} 'Decimal']), 'max')
            Ps.([D{d} 'Decimal']) = findMaxDecLen(get(Axs(j), [D{d} 'TickLabel']));
        end
        resetDecimal(Axs, D{d}, Ps.([D{d} 'Decimal']));
    end
    set(Axs(j), 'Position', PriorPositions{j});
end

%Finds the maximum number of deimals in a TickLabel
function MaxDec = findMaxDecLen(TickLabel)
MaxDec = 0;
for j = 1:length(TickLabel)
    DotIdx = find(TickLabel{j} == '.', 1, 'last');
    if ~isempty(DotIdx)
        DecLen = length(TickLabel{j}) - DotIdx;
        MaxDec = max(DecLen, MaxDec);
    end
end

%Resets the decimal places in the X,Y,Z labels
function resetDecimal(Axs, Direction, Decimal)
D = upper(Direction(1));
if Decimal >= 0
    for j = 1:length(Axs)      
        StrFmt = ['%0.' num2str(Decimal) 'f'];
        StrFunc = str2func([lower(D) 'tickformat']);
        StrFunc(Axs(1), StrFmt);
    end
end
