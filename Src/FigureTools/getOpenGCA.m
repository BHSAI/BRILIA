%getOpenGCA will get a current axes handle ONLY if there is one open.
%Otherwise, will return an empty matrix. Better than simple 'gca' which
%will create a figure and axes if there is none.
%
%  Ax = getOpenGCA
%
%  Ax = getOpenGCA(Gx)
%
%  Ax = getOpenGCA('all')
%
%  [Ax, varargout] = getOpenGCA(varargin{:})
%
%  INPUT
%    Gx: figure or axes handle
%    'all': will get all open axes on all figures
%    varargin: 1xM cell of inputs, in which this function checks if the
%      first one is a graphics handle or 'all'
%
%  OUTPUT
%    Ax: Axes handle(s) that is (are) open or 1st child of figure handle
%    varargout: same as varargin, EXCEPT IF the first varargin cell is a
%      axes handle, then varargout removes the first cell

function [Ax, varargout] = getOpenGCA(varargin)
Ax = [];
varargout{1} = varargin;

if isempty(varargin) 
    if ~isempty(findobj(0, 'type', 'axes'))
        Ax = gca;
    end
    return;
end

if isempty(varargin{1})
    if ~isempty(findobj(0, 'type', 'axes'))
        Ax = gca;
    end
    varargout{1}(1) = [];
    return;
end

if ischar(varargin{1}) && strcmpi(varargin{1}, 'all')
    Ax = findobj(0, 'type', 'axes');
    if isempty(Ax)
        Ax = [];
    end
    varargout{1}(1) = [];
    return;
end    

if strcmpi(class(varargin{1}), 'matlab.graphics.Graphics') %Determine all axes handles
    Ox = varargin{1};
    GetIdx = zeros(length(Ox), 1, 'logical');
    for j = 1:length(Ox)
        if strcmpi(class(Ox(j)), 'matlab.graphics.axis.Axes')
            GetIdx(j) = 1;
        end
    end
    Ax = Ox(GetIdx);
    varargout{1}(1) = [];
    return;
end

if strcmpi(class(varargin{1}(1)), 'matlab.graphics.axis.Axes')
    Ax = varargin{1};
    varargout{1}(1) = [];
    return;
end

if strcmpi(class(varargin{1}(1)), 'matlab.ui.Figure')
    Ax = get(varargin{1}(1), 'children');
    varargout{1}(1) = [];
    return;
end
