%Simple function for showing status of BRILIA if there is a GUI or not.
function showStatus(Msg, varargin)
if isempty(varargin)
    fprintf('Status: %s\n', Msg);
else
    if strcmpi(class(varargin{1}),'matlab.ui.control.UIControl')
        set(varargin{1},'String',Msg);        
        drawnow; %Required short pause to update gui. Otherwise, will appear frozen.
    else
        fprintf('Status: %s\n', Msg);
    end
end