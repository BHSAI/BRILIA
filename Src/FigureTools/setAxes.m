%setAxes will set all the axes fields given the Ax axes handle and the
%param-value inputs. To adjust the title property, add 'Title' in the
%parameter name, (eg, 'TitleFontSize' to adjust title font size, and
%'FontSize' to adjust axes font size.).
%
%  setAxes(Ax, Param, Value, ...)
%
%  setAxes(Gx, Param, Value, ...)
%
%  INPUT
%    Ax: axes handle to modify
%    Gx: figure handle, in which all axes within Gx will be modified
%    Param: any axes property that can be modified by Matlab
%    Value: value to use for the axes property
%
%  OUTPUT
%    Modifies the property of an individual axes or all axes in a figure.
%
%  NOTE
%    To modify title properties, add the string 'Title' before the
%    property, such as 'TitleFontSize'. This will adjust the title axes
%    directly.
%
%  EXAMPLE
%    loadDemoPlot
%    setAxes(gcf, 'XColor', [0 0 0], 'YColor', [0 0 0]);
%    setAxes(gcf, 'TitleString', 'Title', 'TitleFontSize', 14, 'TitleColor', 'r')
%
function setAxes(varargin)
[Axs, varargin] = getOpenGCA(varargin{:});

%Make sure you are given a valid param-value paired system
TitleParamLoc = zeros(length(varargin), 1, 'logical'); %Location of title modifiers
for j = 1:2:length(varargin)
    if ~ischar(varargin{j}) %Not a string param
        error('%s: The varargin is not a perfect param-value pair. Recheck input.', mfilename);
    else
        if ~isempty(min(strfind(lower(varargin{j}), 'title')) == 1)
            varargin{j} = varargin{j}(length('title')+1:end);
            try
                TitleParamLoc(j:j+1) = 1;
            catch
                error('%s: The varargin for a Title parameter is missing its value.', mfilename);
            end
        end
    end
end
TitleVarargin = varargin(TitleParamLoc);
varargin(TitleParamLoc) = [];

for a = 1:length(Axs)
    %Modify the Axes handle
    FieldNames = fieldnames(Axs(a));
    for j = 1:2:length(varargin)
        try
            if startsWith(varargin{j}, {'XLabel', 'YLabel', 'ZLabel'}) %For relabeling the text of axis, ex: XLabel-String, XLabel-FontName, etc
                SpecialOpt = strsplit(varargin{j}, '-');
                if numel(SpecialOpt) == 1
                    SpecialOpt{2} = 'String';
                end
                set(get(Axs(a), SpecialOpt{1}), SpecialOpt{2}, varargin{j+1});                
            elseif any(strcmpi(FieldNames, varargin{j})) %Oddly, isfield does NOT work for axes handles like isfield(gca, 'Units') will fail.
                set(Axs(a), varargin{j}, varargin{j+1});
            end
        catch ME
            warning('%s: Could not set axes properties.\n%s', mfilename, ME.message);
        end
    end
    
    %Modify the Title handle
    if ~isempty(TitleVarargin)
        Tx = get(Axs(a), 'Title');
        for t = 1:2:length(TitleVarargin)
            if strcmpi(TitleVarargin{t}, 'FontSize') %Need to modify the axes TitleFontSizeMultiplier
                AxesFontSize = get(Axs(a), 'FontSize');
                TitleFontSize = TitleVarargin{t + 1};
                Multiplier = TitleFontSize / AxesFontSize;
                set(Axs(a), 'TitleFontSizeMultiplier', Multiplier);
            else
                try
                    set(Tx, TitleVarargin{t}, TitleVarargin{t+1});
                catch
                    warning('%s: Could not set %s for title', mfilename, TitleVarargin{t});
                end
            end
        end
    end
end
