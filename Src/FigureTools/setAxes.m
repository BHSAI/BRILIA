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
StructLoc = cellfun(@isstruct, varargin);
StructVar = cellfun(@(x) struct2array(x, 'input'), varargin(StructLoc), 'un', 0);
varargin = horzcat(varargin(~StructLoc), StructVar{:});

for a = 1:length(Axs)
    %Modify the Axes handle
    FieldNames = fieldnames(Axs(a));
    for j = 1:2:length(varargin)
        try
            if startsWith(varargin{j}, {'XLabel', 'YLabel', 'ZLabel'}) %For relabeling the text of axis, ex: XLabel-String, XLabel-FontName, etc
                SpecialOpt = strsplit(varargin{j}, {'=', '-', ':'});
                if numel(SpecialOpt) == 1
                    SpecialOpt{2} = 'String';
                end
                set(get(Axs(a), SpecialOpt{1}), SpecialOpt{2}, varargin{j+1});
            elseif startsWith(varargin{j}, 'Title') %For relabeling the text of axis, ex: XLabel-String, XLabel-FontName, etc
                SpecialOpt = strsplit(varargin{j}, {'=', '-', ':'});
                if numel(SpecialOpt) == 1
                    SpecialOpt{2} = 'String';
                end
                set(get(Axs(a), SpecialOpt{1}), SpecialOpt{2}, varargin{j+1});
            elseif any(strcmpi(FieldNames, varargin{j})) %CODING_NOTE: Using strcmpmmi instead of isfield because oddly, isfield does NOT work for axes handles like isfield(gca, 'Units') will fail.
                set(Axs(a), varargin{j}, varargin{j+1});
            end
        catch ME
            error('%s: Could not set axes properties.\n%s', mfilename, ME.message);
        end
    end
end
