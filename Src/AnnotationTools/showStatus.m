%Simple function for showing status of BRILIA if there is a GUI or not.
function showStatus(Msg,varargin)
Msg = strrep(Msg, '\', '\\');
if isempty(varargin)
    fprintf([Msg '\n']);
else
    if strcmpi(class(varargin{1}),'matlab.ui.control.UIControl')
        set(varargin{1},'String',Msg);        
        pause(0.01); %Required short pause to update gui. Otherwise, will appear frozen.
    else
        fprintf([Msg '\n']);
    end
end
