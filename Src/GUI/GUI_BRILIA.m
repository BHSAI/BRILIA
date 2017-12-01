function varargout = GUI_BRILIA(varargin)
% GUI_BRILIA MATLAB code for GUI_BRILIA.fig
%      GUI_BRILIA, by itself, creates a new GUI_BRILIA or raises the existing
%      singleton*.
%
%      H = GUI_BRILIA returns the handle to a new GUI_BRILIA or the handle to
%      the existing singleton*.
%
%      GUI_BRILIA('CALLBACK', hObject, eventData, handles, ...) calls the local
%      function named CALLBACK in GUI_BRILIA.M with the given input arguments.
%
%      GUI_BRILIA('Property', 'Value', ...) creates a new GUI_BRILIA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_BRILIA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_BRILIA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_BRILIA

% Last Modified by GUIDE v2.5 24-Nov-2017 17:04:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',      mfilename, ...
                   'gui_Singleton', gui_Singleton, ...
                   'gui_OpeningFcn',@GUI_BRILIA_OpeningFcn, ...
                   'gui_OutputFcn', @GUI_BRILIA_OutputFcn, ...
                   'gui_LayoutFcn', [], ...
                   'gui_Callback',  []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function GUI_BRILIA_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.text_Version, 'String', ['Version ' BRILIA('getversion')]);

SpeciesList = getGeneDatabase('getlist', 'suppress');
SortOrder = {'mouse', 'human'};
for k = 1:length(SortOrder)
    Loc = find(contains(SpeciesList, SortOrder{k}, 'IgnoreCase', true));
    SpeciesList([k Loc]) = SpeciesList([Loc k]);
end
set(handles.popupmenu_Species, 'String', SpeciesList, 'Value', 1);
popupmenu_Species_Callback(hObject, eventdata, handles)

handles.SaveFileNames = {}; %List of output files
handles.output = hObject;
guidata(hObject, handles);

function varargout = GUI_BRILIA_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
%All popup lists

function popupmenu_Species_Callback(hObject, eventdata, handles)
%Need to ensure the the strain list and selection matches with species
SpeciesList = get(handles.popupmenu_Species, 'string');
SpeciesLoc  = get(handles.popupmenu_Species, 'value');
StrainList  = getStrainList(SpeciesList{SpeciesLoc});
StrainNum   = get(handles.popupmenu_Strain, 'value');
if StrainNum > length(StrainList)
    set(handles.popupmenu_Strain, 'value', 1);
end
set(handles.popupmenu_Strain, 'string', StrainList); 
function popupmenu_Species_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function popupmenu_Strain_Callback(hObject, eventdata, handles)
function popupmenu_Strain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function popupmenu_Chain_Callback(hObject, eventdata, handles)
function popupmenu_Chain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function popupmenu_FileType_Callback(hObject, eventdata, handles)
function popupmenu_FileType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'), get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

%--------------------------------------------------------------------------
%All checkboxes

function checkbox_Vfunctional_Callback(hObject, eventdata, handles)

function checkbox_Vpseudo_Callback(hObject, eventdata, handles)

function checkbox_Vorf_Callback(hObject, eventdata, handles)

function checkbox_Dforward_Callback(hObject, eventdata, handles)

function checkbox_Dinverse_Callback(hObject, eventdata, handles)

function checkbox_CheckSeqDirection_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
%All editable text

function edit_DevPerc_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 0, 100, 0, 'double');
if Status == 0
    set(hObject, 'String', ValidStr);
end
function edit_DevPerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'),  get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

function edit_NumProc_Callback(hObject, eventdata, handles)
CheckStr = get(hObject, 'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr, 1, feature('numCores'), 1, 'int');
if Status == 0
    set(hObject, 'String', ValidStr);
end
function edit_NumProc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject, 'BackgroundColor'),  get(0, 'defaultUicontrolBackgroundColor'))
    set(hObject, 'BackgroundColor', 'white');
end

%--------------------------------------------------------------------------
%All buttons

function pushbutton_OpenSettingFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uigetfile('*.txt', 'Open the setting text file for BRILIA');
if isnumeric(FilePath); return; end
set(handles.text_SettingFile, 'String', [FilePath FileName]);

P = BRILIA('getinput'); %Return the current defaults in BRILIA. Use as input to readSettingFile, to ensure P is the update with BRILIA version.
P = readSettingFile([FilePath FileName], P); %Load here since you want it to update all the fields in GUI too.

%Update the Species
SpeciesList = get(handles.popupmenu_Species, 'string');
SpeciesLoc = findCell(SpeciesList, P.Species, 'MatchCase', 'any');
if SpeciesLoc(1) > 0 %Update the Species and Strain list
    set(handles.popupmenu_Species, 'value', SpeciesLoc); %Update the selection for species   
    StrainList = getStrainList(SpeciesList{SpeciesLoc}); %Update the strain list for this species
    set(handles.popupmenu_Strain, 'string', StrainList);
end

%Update the Strain
SpeciesNum  = get(handles.popupmenu_Species, 'value');
SpeciesName = SpeciesList{SpeciesNum};
StrainList  = getStrainList(SpeciesName);
if ~isempty(P.Strain)
    StrainNum = findCell(StrainList, P.Strain, 'MatchCase', 'any', 'MatchWord', 'partial');
    if StrainNum(1) == 0
        StrainNum = 1;
    end
else
    StrainNum = 1; %all
end
set(handles.popupmenu_Strain, 'value', StrainNum);

%Update the clustering deviation %
if ~isempty(P.DevPerc)
    set(handles.edit_DevPerc, 'string', num2str(P.DevPerc));
    edit_DevPerc_Callback(handles.edit_DevPerc, eventdata, handles); %Recall to ensure valid value
end

%Update the checkbox for V gene functions
if ~isempty(P.Vfunction)
    if strcmpi(P.Vfunction, 'all')
        VfunctCell = {'f', 'p', 'o'}; %do all of it.
    else
        VfunctCell = regexpi(strrep(P.Vfunction, ' ', ''), ',', 'split');
    end
    NewVcheckbox = zeros(1, 3); %See if setting file is valid
    if ~isempty(VfunctCell) %Just make sure regexp doesn't yield empty
        for k = 1:length(VfunctCell)
            switch lower(VfunctCell{k}(1))
                case 'f' %Functional
                    NewVcheckbox(1) = 1;
                case 'p' %Pseudogene
                    NewVcheckbox(2) = 1;
                case 'o' %Open reading frame
                    NewVcheckbox(3) = 1;
            end
        end
    end
    
    %Make sure at least 1 check box is selected, then update
    if max(NewVcheckbox) > 0
        set(handles.checkbox_Vfunctional, 'value', NewVcheckbox(1));
        set(handles.checkbox_Vpseudo, 'value', NewVcheckbox(2));
        set(handles.checkbox_Vorf, 'value', NewVcheckbox(3));                
    end
end

%Update the checkbox for the D gene directions
if ~isempty(P.Ddirection)
    if strcmpi(P.Ddirection, 'all')
        DdirCell = {'f', 'i'}; %do forward and inverse
    else
        DdirCell = regexpi(strrep(P.Ddirection, ' ', ''), ',', 'split');
    end
    NewDcheckbox = zeros(1, 2);
    for k = 1:length(DdirCell)
        switch lower(DdirCell{k}(1))
            case 'f' %Functional
                NewDcheckbox(1) = 1;
            case 'i' %Pseudogene
                NewDcheckbox(2) = 1;
        end
    end
    
    %Make sure at least 1 check box is selected, then update
    if max(NewVcheckbox) > 0
        set(handles.checkbox_Dforward, 'value', NewDcheckbox(1));
        set(handles.checkbox_Dinverse, 'value', NewDcheckbox(2));
    end
end

function pushbutton_InputFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uigetfile('*.fa*;*.*sv', 'Open the sequence file to process.');
if isnumeric(FilePath); return; end
InputFileName = fullfile(FilePath, FileName);
set(handles.text_InputFileName, 'string', InputFileName);

%Generate a default output file name
[OutputFilePath, OutputFileName, ~] = parseFileName(InputFileName);
DotLoc = find(OutputFileName == '.');
Version = BRILIA('getversion');
SaveFilePath = fullfile(OutputFilePath, OutputFileName(1:DotLoc(end)-1));
SaveFileName = [OutputFileName(1:DotLoc(end)-1) '.BRILIAv' Version(1) '.csv'];
OutputFileName = fullfile(SaveFilePath, SaveFileName); 
set(handles.text_OutputFileName, 'string', OutputFileName);

function pushbutton_OutputFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uiputfile('*.fa*;*.xls*;*.csv', 'Save the results to here.');
if isnumeric(FilePath); return; end
OutputFileName = fullfile(FilePath, FileName);
set(handles.text_OutputFileName, 'string', OutputFileName);

function pushbutton_Start_Callback(hObject, eventdata, handles)
P = BRILIA('getinput');
%--------------------------------------------------------------------------
%Ensure at least 1 checkmark for Vfunction, then set it
Vbox = zeros(1, 3, 'logical');
Vbox(1) = get(handles.checkbox_Vfunctional, 'value');
Vbox(2) = get(handles.checkbox_Vpseudo, 'value');
Vbox(3) = get(handles.checkbox_Vorf, 'value');
if max(Vbox) == 0 
    Msg = 'Must select 1 V function checkbox.';
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 0 0]);
    return
else 
    %Format V function input for BRILIA
    VboxStr = {'f', 'p', 'orf'}; %BRILIA-recognized string inputs, in order of Vbox
    if sum(Vbox) == length(Vbox)
        P.Vfunction = 'all';
    else
        VboxStr = VboxStr(Vbox);
        P.Vfunction = sprintf(repmat('%s, ', 1, length(VboxStr)), VboxStr{:});
        P.Vfunction(end) = [];
    end

    Msg = '';
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 1 1]);
end

%--------------------------------------------------------------------------
%Ensure at least 1 checkmark for Ddirection, then set it
Dbox = zeros(1, 2, 'logical');
Dbox(1) = get(handles.checkbox_Dforward, 'value');
Dbox(2) = get(handles.checkbox_Dinverse, 'value');
if max(Dbox) == 0 
    Msg = 'Must select 1 D direction checkbox.';
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 0 0]);
    return
else 
    %Format D direction input for BRILIA
    DboxStr = {'fwd', 'inv'}; %BRILIA-recognized string inputs, in order of Dbox
    if sum(Dbox) == length(Dbox)
        P.Ddirection = 'all';
    else
        P.Ddirection = DboxStr{Dbox > 0};
    end
    
    Msg = '';
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 1 1]);    
end

CheckDirBox = get(handles.checkbox_CheckSeqDirection, 'value');
if CheckDirBox
    P.CheckSeqDir = 'y';
else
    P.CheckSeqDir = 'n';
end
%--------------------------------------------------------------------------
%Ensure file name exists, then set it
FullFileName = get(handles.text_InputFileName, 'String');
if isempty(FullFileName)
    Msg = 'Select the sequence file.';
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 0 0]);    
    return
else
    Msg = sprintf('Opening [ %s ]',  FullFileName);
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 1 1]);        
end
if ~exist(FullFileName, 'file')
    Msg = sprintf('Input file does not exists at [ %s ].', FullFileName);
    set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 0 0]);    
    return
end
P.InputFile = FullFileName;

P.OutputFile = get(handles.text_OutputFileName, 'String');

%Extract other information
P.DevPerc = str2double(get(handles.edit_DevPerc, 'String'));
P.StatusHandle = handles.text_Status;

StrainList = get(handles.popupmenu_Strain, 'String');
P.Strain = StrainList{get(handles.popupmenu_Strain, 'Value')};

SpeciesList = get(handles.popupmenu_Species, 'String');
P.Species = SpeciesList{get(handles.popupmenu_Species, 'Value')};

%--------------------------------------------------------------------------
%Open up the processors now
Msg = 'Opening up processors for BRILIA';
set(handles.text_Status, 'String', Msg, 'ForegroundColor', [1 1 1]);        
drawnow;
try 
    P.NumProc = convStr2Num(get(handles.edit_NumProc, 'String'));
catch
    P.NumProc = 1;
end
setCores(P.NumProc);
Msg = sprintf('Using %d Processors', getActiveCores);
set(handles.text_ParallelProc, 'String', Msg);
drawnow; %Ensure text updates

%Begin BRILIA
Msg = 'Starting BRILIA';
set(handles.text_Status, 'String', Msg, 'ForegroundColor', [0 0.8 0]);
drawnow;

SaveFileNames = BRILIA(FullFileName, 'SuppressIntro', 'y', P);
RunTime = toc;

Msg = sprintf('Job completed in %1.2f min.\n', RunTime/60);
set(handles.text_Status, 'String', Msg, 'ForegroundColor', [0 0.8 0]);

handles.SaveFileNames = SaveFileNames;
guidata(hObject, handles);

function pushbutton_PlotTree_Callback(hObject, eventdata, handles)
if ~isempty(handles.SaveFileNames)
    SaveFileName = handles.SaveFileNames{1}; %Only get the first one
    GUI_plotTree([], SaveFileName);
else
    GUI_plotTree;    
end

%--------------------------------------------------------------------------
%Special functions for the GUI

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

%This store the list of strain names
function StrainList = getStrainList(SpeciesName)
DB = getGeneDatabase(SpeciesName, 'suppress');
StrainList = getUnqStrain(DB, 4);
if ~strcmpi(StrainList{1}, 'All')
    StrainList = cat(1, 'All', StrainList);
end
