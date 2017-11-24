function varargout = GUI_plotTree(varargin)
% GUI_PLOTTREE MATLAB code for GUI_plotTree.fig
%      GUI_PLOTTREE, by itself, creates a new GUI_PLOTTREE or raises the existing
%      singleton*.
%
%      H = GUI_PLOTTREE returns the handle to a new GUI_PLOTTREE or the handle to
%      the existing singleton*.
%
%      GUI_PLOTTREE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PLOTTREE.M with the given input arguments.
%
%      GUI_PLOTTREE('Property', 'Value', ...) creates a new GUI_PLOTTREE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_plotTree_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_plotTree_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_plotTree

% Last Modified by GUIDE v2.5 22-Aug-2017 10:37:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_plotTree_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_plotTree_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function GUI_plotTree_OpeningFcn(hObject, eventdata, handles, varargin)
Failed = false; 
if length(varargin) >= 2 && ischar(varargin{2}) 
    try 
        handles.FileName = varargin{2};
        [handles.VDJdata, handles.VDJheader] = openSeqData(handles.FileName);
        set(handles.text_FileName, 'String', handles.FileName);
    catch
        Failed = true;
    end
end
if Failed
    handles.FileName = [];
    handles.VDJdata  = {};
    handles.VDJheader = {};
end    
handles.output = hObject;
guidata(hObject, handles);

function varargout = GUI_plotTree_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pushbutton_OpenFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uigetfile('*.csv*;*.xls*', 'Open the processed sequence file from BRILIA');
if isnumeric(FilePath); return; end
handles.FileName = fullfile(FilePath, FileName);
try
    [handles.VDJdata, handles.VDJheader] = openSeqData(handles.FileName);
catch
    return
end
set(handles.text_FileName, 'string', handles.FileName);
guidata(hObject, handles);

function pushbutton_Plot_Callback(hObject, eventdata, handles)
[Gx, ~] = getOpenPlotHandles;
if ~isempty(Gx)
    close(Gx{:});
end

try 
    P = plotTree('getinput');
    
    P.GetGrpNum = str2double(get(handles.edit_GetGrpNum, 'String'));
    if P.GetGrpNum == 0; P.GetGrpNum = []; end
    
    P.GetSeqNum = str2double(get(handles.edit_GetSeqNum, 'String'));
    if P.GetSeqNum == 0; P.GetSeqNum = []; end

    P.GetSizeRange = [str2double(get(handles.edit_GetSizeRange1, 'String')) ...
                      str2double(get(handles.edit_GetSizeRange2, 'String')) ];
    if P.GetSizeRange(end) == 0; P.GetSizeRange = []; end
        
    P.GetCDR3seq = get(handles.edit_GetCDR3seq, 'String');
    
    P.ShowOnly = true;
    
    plotTree(handles.VDJdata, handles.VDJheader, P)
catch
end
    
function edit_GetGrpNum_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 0, inf, 0, 'int');
if Status == 0
    set(hObject, 'String', ValidStr);
end
function edit_GetGrpNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function edit_GetSeqNum_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 0, inf, 0, 'int');
if Status == 0
    set(hObject, 'String', ValidStr);
end
function edit_GetSeqNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function edit_GetSizeRange1_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 0, inf, 0, 'int');
if Status == 0
    set(hObject, 'String', ValidStr);
end

%Ensure that GetSizeRange2 is not bigger
Value1 = str2double(get(handles.edit_GetSizeRange1, 'String'));
Value2 = str2double(get(handles.edit_GetSizeRange2, 'String'));
if Value1 > Value2 %Need to reset the OTHER text field
    set(handles.edit_GetSizeRange2, 'String', num2str(Value1));
end
function edit_GetSizeRange1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function edit_GetSizeRange2_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 0, inf, 0, 'int');
if Status == 0
    set(hObject, 'String', ValidStr);
end

%Ensure that GetSizeRange1 is not bigger
Value1 = str2double(get(handles.edit_GetSizeRange1, 'String'));
Value2 = str2double(get(handles.edit_GetSizeRange2, 'String'));
if Value1 > Value2 %Need to reset the OTHER text field
    set(handles.edit_GetSizeRange1, 'String', num2str(Value2));
end
function edit_GetSizeRange2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function edit_GetCDR3seq_Callback(hObject, eventdata, handles)
CurStr = upper(get(hObject, 'String'));
BadPattern = ['[^' int2aa(1:20) 'X]'];
BadCharLoc = regexp(CurStr, BadPattern);
CurStr(BadCharLoc) = 'X';
set(hObject, 'String', CurStr);
function edit_GetCDR3seq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

%Checks to see if string is a valid number and within a range
function [Status, ValidStr] = checkStringIsNumeric(CheckStr, Min, Max, Default, varargin)
Type = 'double';
if ~isempty(varargin)
    Type = varargin{1};
end

%Start by assuming it is valid
ValidStr = CheckStr;
Status = 1;

%Make sure there are up to 1 decimal
DotLoc = find(CheckStr == '.');
if length(DotLoc) > 1 %Invalid, too many dots.
    ValidStr = num2str(Default);
    Status = 0;
    return
end

%Make sure there are no non-digit values
if ~isempty(regexpi(CheckStr, '[^\d\.]', 'once'))
    ValidStr = num2str(Default);
    Status = 0;
    return
end

%Make sure the value type is correct
CheckNum = str2double(CheckStr);
if strcmpi(Type, 'int')
    if mod(CheckNum, 1) ~= 0 %Has a decimal
        CheckNum = round(CheckNum);
        Status = 0;
    end
end

%Make sure the value range is correct
if CheckNum < Min
    CheckNum = Min;
    Status = 0;
elseif CheckNum > Max
    CheckNum = Max;
    Status = 0;
end

%Return the valid number as a string
if Status == 0
    ValidStr = num2str(CheckNum);
end

function [Gx, Ax] = getOpenPlotHandles
OpenFigures = findobj('Type', 'figure');
GuiFigLoc = zeros(length(OpenFigures), 1, 'logical');
for j = 1:length(OpenFigures)
    if contains(lower(OpenFigures(j).Name), 'gui') %file name can't have gui in there.
        GuiFigLoc(j) = 1;
    end
end
OpenFigures(GuiFigLoc) = [];

%Return a cell array of figure and axes handles
Gx = cell(size(OpenFigures));
Ax = cell(size(OpenFigures));
for j = 1:length(OpenFigures)
    Gx{j} = OpenFigures(j);
    Ax{j} = get(OpenFigures(j), 'CurrentAxes');
end

function pushbutton_SavePlots_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uiputfile('*.tif*;*.png*;*.jpg;', 'Save images as');
if isnumeric(FilePath); return; end

[SavePath, SaveFile, SaveExt] = parseFileName(fullfile(FilePath, FileName));
if isempty(SaveExt) %Then set default to tif
    SaveExt = '.tif';
end

DotLoc = find(SaveFile == '.');
if isempty(DotLoc)
    SaveName = fullfile(SavePath, [SaveFile SaveExt]);
else
    SaveName = fullfile([SaveFile(1:DotLoc-1) SaveExt]);
end

%Get the title of the plot, which has the tree name
[Gx, Ax] = getOpenPlotHandles;
for j = 1:length(Gx)
    CurGx = Gx{j};
    CurAx = Ax{j};
    
    %Attempt to get the grp number info
    TitleStr = get(get(CurAx, 'Title'), 'String');
    CurGrpStr = regexp([TitleStr{:}], 'Grp\s+(\d+),', 'tokens');
    if ~isempty(CurGrpStr)
        CurGrpStr = ['Grp' CurGrpStr{1}{1}];
    else
        CurGrpStr = ['A' num2str(j)]; %Assign something
    end
    
    savePlot(CurGx, 'DPI', 300, 'SaveAs', SaveName, 'SavePre', CurGrpStr);
    
    %Saving the image
    %sprintf('%s%s.Grp%s%s', SavePath, SaveNamePre, CurGrpStr, SaveExt);
%     DPI = 300;
%     switch lower(SaveExt)
%         case '.tif'
%             print(CurGx, SaveName, '-dtiff', ['-r' num2str(DPI)]);
%         case '.jpg'
%             print(CurGx, SaveName, '-djpeg', ['-r' num2str(DPI)]);
%         case '.png'
%             print(CurGx, SaveName, '-dpng', ['-r' num2str(DPI)]);                
%         case '.fig'
%             saveas(CurGx, SaveName);
%     end
end

function pushbutton_ClosePlots_Callback(hObject, eventdata, handles)
%Close all open figures, and get the others
[Gx, ~] = getOpenPlotHandles;
if ~isempty(Gx)
    close(Gx{:});
end

function figure1_CloseRequestFcn(hObject, eventdata, handles)
%Close all open figures, and get the others
[Gx, ~] = getOpenPlotHandles;
if ~isempty(Gx)
    close(Gx{:});
end
delete(hObject);
