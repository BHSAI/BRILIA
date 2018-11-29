%showStatus will pass on a status message to stdout or to a UI message box.
%
%  showStatus(Msg, StatusHandle)
%
%  INPUT
%    Msg: message string to show 
%    StatusHandle: handle to a UI object with "String" property to set. If 
%      empty, will print message to stdout.
%
function showStatus(Msg, varargin)
if isempty(varargin)
    fprintf('Status: %s\n', Msg);
else
    if strcmpi(class(varargin{1}),'matlab.ui.control.UIControl')
        set(varargin{1}, 'String', Msg);        
        drawnow('nocallbacks');
    else
        fprintf('Status: %s\n', Msg);
    end
end