function varargout = BRILIAgui(varargin)
% BRILIAGUI MATLAB code for BRILIAgui.fig
%      BRILIAGUI, by itself, creates a new BRILIAGUI or raises the existing
%      singleton*.
%
%      H = BRILIAGUI returns the handle to a new BRILIAGUI or the handle to
%      the existing singleton*.
%
%      BRILIAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRILIAGUI.M with the given input arguments.
%
%      BRILIAGUI('Property','Value',...) creates a new BRILIAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BRILIAgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BRILIAgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BRILIAgui

% Last Modified by GUIDE v2.5 08-Feb-2017 14:19:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BRILIAgui_OpeningFcn, ...
                   'gui_OutputFcn',  @BRILIAgui_OutputFcn, ...
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

function BRILIAgui_OpeningFcn(hObject, eventdata, handles, varargin)
%Make sure BRILIA paths are added correctly
CurPaths = regexp(path,';','split')';
MainPath = mfilename('fullpath');
SlashLoc = regexp(MainPath,'\\|\/');
MainPath = MainPath(1:SlashLoc(end-1)-1); %Note that matlab does not save the final slash. you do SlashLoc(end-1) to get to the main BRILIA folder
HavePath = 0;
for p = 1:length(CurPaths)
    if strcmp(CurPaths{p},MainPath(1:end)) 
        HavePath = 1;
        break
    end
end
if HavePath == 0 %Matlab doesn't have path, so must add it
    disp('Adding BRILIA paths to MatLab');
    addpath(genpath(MainPath));
end

handles.SaveFileNames = {}; %List of output files

handles.output = hObject;
guidata(hObject, handles);

function varargout = BRILIAgui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%--------------------------------------------------------------------------
%All popup lists

function popupmenu_Species_Callback(hObject, eventdata, handles)
%Need to ensure the the strain list and selection matches with species
SpeciesList = get(handles.popupmenu_Species,'string');
SpeciesLoc = get(handles.popupmenu_Species,'value');
StrainList = getStrainList(SpeciesList{SpeciesLoc});
StrainNum = get(handles.popupmenu_Strain,'value');
if StrainNum > length(StrainList)
    set(handles.popupmenu_Strain,'value',1);
end
set(handles.popupmenu_Strain,'string',StrainList); 
function popupmenu_Species_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_Strain_Callback(hObject, eventdata, handles)
function popupmenu_Strain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_Chain_Callback(hObject, eventdata, handles)
function popupmenu_Chain_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_FileType_Callback(hObject, eventdata, handles)
function popupmenu_FileType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%All checkboxes

function checkbox_Vfunctional_Callback(hObject, eventdata, handles)

function checkbox_Vpseudo_Callback(hObject, eventdata, handles)

function checkbox_Vorf_Callback(hObject, eventdata, handles)

function checkbox_Dforward_Callback(hObject, eventdata, handles)

function checkbox_Dinverse_Callback(hObject, eventdata, handles)

%--------------------------------------------------------------------------
%All editable text

function edit_DevPerc_Callback(hObject, eventdata, handles)
CheckStr = get(hObject,'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr,0,100,0,'double');
if Status == 0
    set(hObject,'String',ValidStr);
end
function edit_DevPerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_NumProc_Callback(hObject, eventdata, handles)
CheckStr = get(hObject,'String');
[Status, ValidStr] = checkStringIsNumeric(CheckStr,1,feature('numCores'),1,'int');
if Status == 0
    set(hObject,'String',ValidStr);
end
function edit_NumProc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%All buttons

function pushbutton_OpenSettingFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uigetfile('*.txt','Open the setting text file for BRILIA');
if isnumeric(FilePath) && FilePath == 0 %User did not open file. return.
    return
end
set(handles.text_SettingFile,'String',[FilePath FileName]);

P = BRILIA('getinput'); %Return the current defaults in BRILIA. Use as input to readSettingFile, to ensure P is the update with BRILIA version.
P = readSettingFile([FilePath FileName],P); %Load here since you want it to update all the fields in GUI too.

%Update the GUI settings now

%Update the Species selection
SpeciesList = get(handles.popupmenu_Species,'string');
SpeciesLoc = findCell(SpeciesList,P.Species,'MatchCase','any');
if SpeciesLoc(1) > 0 %Update the Species and Strain list
    set(handles.popupmenu_Species,'value',SpeciesLoc); %Update the selection for species   
    StrainList = getStrainList(SpeciesList{SpeciesLoc}); %Update the strain list for this species
    set(handles.popupmenu_Strain,'string',StrainList);
end

%Update the strain selection
SpeciesNum = get(handles.popupmenu_Species,'value');
SpeciesName = SpeciesList{SpeciesNum};
StrainList = getStrainList(SpeciesName);
if ~isempty(P.Strain)
    StrainNum = findCell(StrainList,P.Strain,'MatchCase','any','MatchWord','partial');
    if StrainNum(1) == 0
        StrainNum = 1;
    end
else
    StrainNum = 1; %all
end
set(handles.popupmenu_Strain,'value',StrainNum);

%Update the clustering deviation %
if ~isempty(P.DevPerc)
    set(handles.edit_DevPerc,'string',num2str(P.DevPerc));
    edit_DevPerc_Callback(handles.edit_DevPerc, eventdata, handles); %Recall to ensure valid value
end

%Update the checkbox for V gene functions
if ~isempty(P.Vfunction)
    if strcmpi(P.Vfunction,'all')
        VfunctCell = {'f','p','o'}; %do all of ti.
    else
        VfunctCell = regexpi(P.Vfunction,',','split');
    end
    NewVcheckbox = zeros(1,3); %See if setting file is valid
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
        set(handles.checkbox_Vfunctional,'value',NewVcheckbox(1));
        set(handles.checkbox_Vpseudo,'value',NewVcheckbox(2));
        set(handles.checkbox_Vorf,'value',NewVcheckbox(3));                
    end
end

%Update the checkbox for the D gene directions
if ~isempty(P.Ddirection)
    if strcmpi(P.Ddirection,'all')
        DdirCell = {'f','i'}; %do forward and inverse
    else
        DdirCell = regexpi(P.Ddirection,',','split');
    end
    NewDcheckbox = zeros(1,2);
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
        set(handles.checkbox_Dforward,'value',NewDcheckbox(1));
        set(handles.checkbox_Dinverse,'value',NewDcheckbox(2));
    end
end

function pushbutton_OpenFile_Callback(hObject, eventdata, handles)
[FileName, FilePath] = uigetfile('*.fa*;*.xls*;*.csv','Open the sequence file to process');
if isnumeric(FilePath) && FilePath == 0 %User did not open file. return.
    return
end

%Determine the file type
DotLoc = find(FileName == '.');
FileExt = FileName(DotLoc(end):end);
if strcmpi(FileExt,'.fa') || strcmpi(FileExt,'.fasta')
    FileType = 'fasta';
elseif strcmpi(FileExt,'.fastq')
    FileType = 'fastq';
elseif strcmpi(FileExt(1:4),'.xls')
    FileType = 'excel';
else %Default to delimited if undetermined or if it really is it.
    FileType = 'delimited';
end

%Update the GUI text box for file name
set(handles.text_FileName,'string',[FilePath FileName]);

%Update the GUI textbox and listbox for files
FileTypeList = get(handles.popupmenu_FileType,'string');
FileTypeLoc = findCell(FileTypeList,FileType,'MatchCase','any','MatchWord','partial');
if FileTypeLoc(1) > 0
    set(handles.popupmenu_FileType,'value',FileTypeLoc(1));
end

function pushbutton_Start_Callback(hObject, eventdata, handles)
P = BRILIA('getinput');
%--------------------------------------------------------------------------
%Ensure at least 1 checkmark for Vfunction, then set it
Vbox = zeros(1,3,'logical');
Vbox(1) = get(handles.checkbox_Vfunctional,'value');
Vbox(2) = get(handles.checkbox_Vpseudo,'value');
Vbox(3) = get(handles.checkbox_Vorf,'value');
if max(Vbox) == 0 
    Msg = 'Must select 1 V function checkbox.';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 0 0]);
    return
else 
    %Format V function input for BRILIA
    VboxStr = {'f','p','orf'}; %BRILIA-recognized string inputs, in order of Vbox
    if sum(Vbox) == length(Vbox);
        P.Vfunction = 'all';
    else
        VboxStr = VboxStr(Vbox);
        P.Vfunction = sprintf(repmat('%s,',1,length(VboxStr)),VboxStr{:});
        P.Vfunction(end) = [];
    end

    Msg = '';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 1 1]);
end

%--------------------------------------------------------------------------
%Ensure at least 1 checkmark for Ddirection, then set it
Dbox = zeros(1,2,'logical');
Dbox(1) = get(handles.checkbox_Dforward,'value');
Dbox(2) = get(handles.checkbox_Dinverse,'value');
if max(Dbox) == 0 
    Msg = 'Must select 1 D direction checkbox.';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 0 0]);
    return
else 
    %Format D direction input for BRILIA
    DboxStr = {'fwd','inv'}; %BRILIA-recognized string inputs, in order of Dbox
    if sum(Dbox) == length(Dbox)
        P.Ddirection = 'all';
    else
        P.Ddirection = DboxStr{Dbox > 0};
    end
    
    Msg = '';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 1 1]);    
end

%--------------------------------------------------------------------------
%Ensure file name exists, then set it
FullFileName = get(handles.text_FileName,'String');
if isempty(FullFileName)
    Msg = 'Select the sequence file.';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 0 0]);    
    return
else
    P.FullFileNames = FullFileName;
    
    Msg = '';
    set(handles.text_Status,'String',Msg,'ForegroundColor',[1 1 1]);        
end

%--------------------------------------------------------------------------
%Determine file type and file delimited
FileTypeList = get(handles.popupmenu_FileType,'String');
FileType = FileTypeList{get(handles.popupmenu_FileType,'value')};
switch FileType
    case 'delimited, semicolon'
        P.Delimiter = ';';
        P.FileType = 'delimited';
    case 'delimited, comma'
        P.Delimiter = ',';
        P.FileType = 'delimited';
    case 'delimited, tab'
        P.Delimiter = '\t';
        P.FileType = 'delimited';
end

%--------------------------------------------------------------------------
%Extract other information
P.DevPerc = str2double(get(handles.edit_DevPerc,'String'));
P.StatusHandle = handles.text_Status;

StrainList = get(handles.popupmenu_Strain,'String');
P.Strain = StrainList{get(handles.popupmenu_Strain,'Value')};

SpeciesList = get(handles.popupmenu_Species,'String');
P.Species = SpeciesList{get(handles.popupmenu_Species,'Value')};

%--------------------------------------------------------------------------
%Open up the processors now
Msg = 'Opening up processors for BRILIA';
set(handles.text_Status,'String',Msg,'ForegroundColor',[1 1 1]);        
pause(0.01);

NumProc = str2double(get(handles.edit_NumProc,'String'));
ps = parallel.Settings;
ps.Pool.AutoCreate = false; %Ensure that parfor is not automatically run.
PoolName = gcp('nocreate');
if isempty(PoolName)
    if NumProc > 1
        parpool(NumProc);
        PoolName = gcp('nocreate');
    end
else
    if PoolName.NumWorkers ~= NumProc
        delete(gcp('nocreate'));
        if NumProc > 1
            parpool(NumProc);
        end
        PoolName = gcp('nocreate');
    end
end

if ~isempty(PoolName) && PoolName.NumWorkers > 1
    Msg = sprintf('Using %d Processors',PoolName.NumWorkers);
    set(handles.text_ParallelProc,'String',Msg);
else
    Msg = sprintf('Using %d Processors',1);
    set(handles.text_ParallelProc,'String',Msg);
end
pause(0.01); %Ensure text updates

%--------------------------------------------------------------------------
%Begin BRILIA
Msg = 'Starting BRILIA. Normally takes around 12 ms per 400 bp seq using 4 processors. Smaller sequences are faster.';
set(handles.text_Status,'String',Msg,'ForegroundColor',[0 0.8 0]);
pause(0.01);

[RunTime, SeqCount, SaveFileNames] = BRILIA(P);

%Finished BRILIA
Msg = sprintf('Job completed.\n  Run Time: %1.2f min.\n  Num of Seq: %d.',RunTime/60, SeqCount);
set(handles.text_Status,'String',Msg,'ForegroundColor',[0 0.8 0]);

handles.SaveFileNames = SaveFileNames;
guidata(hObject, handles);

function pushbutton_PlotTree_Callback(hObject, eventdata, handles)
if ~isempty(handles.SaveFileNames)
    SaveFileName = handles.SaveFileNames{1}; %Only get the first one
    plotTreeDataGUI([],SaveFileName);
else
    plotTreeDataGUI;    
end

%--------------------------------------------------------------------------
%Special functions for the GUI

%Checks to see if string is a valid number and within a range
function [Status, ValidStr] = checkStringIsNumeric(CheckStr,Min,Max,Default,varargin)
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
if ~isempty(regexpi(CheckStr,'[^\d\.]','once'))
    ValidStr = num2str(Default);
    Status = 0;
    return
end

%Make sure the value type is correct
CheckNum = str2double(CheckStr);
if strcmpi(Type,'int')
    if mod(CheckNum,1) ~= 0 %Has a decimal
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
if strcmpi(SpeciesName,'Human')
    StrainList = {'All'};
elseif strcmpi(SpeciesName,'Mouse')
%     %This is how you would get the strain list generally
%     [Vmap,Dmap,Jmap] = getCurrentDatabase('change',SpeciesName);
%     StrainList = [{'All'}; getUnqStrain(Vmap,Dmap,Jmap,'Condense',4)];
%     StrainListText = sprintf(repmat('%s;',1,length(StrainList)),StrainList{:})
%     StrainListText(end) = [];
    StrainListText = 'All;129/Sv;A/J;AKR;BALB.K,BALB/b,BALB/c;C57BL,C57BL/10,C57BL/6,C57BL/6J;C58;CB.20;CBA;CE;CLA-2/Cn;Cloned;DBA/2;I/St;MRL/Mp-lpr,MRL/lpr;NFS;NZB;NZW;RF;RIII;SJL';
    StrainList = regexp(StrainListText,';','split')';
end
