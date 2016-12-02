%buildTree2 will go through VDJdata and plot the trees per unique group.
%HOWEVER, but will first look for all unique sequences in both files,
%assemble the color code scheme, and then apply them throughout the entire
%trees for all files. This is useful if you have a main tree, subtrees from
%2 files, and you want to color code the CDR3 regions the same way.

%1) Identify lowest mutated V and J and D.
%2) Assemble reference sequence
%3) Calculated distances from ref to seq
%4) Calculated distances from seq to seq1
%etc etc

function VDJdata = buildTree2(varargin)
P = inputParser;
addOptional(P,'VDJdata',[],@(x) iscell(x) || ischar(x));
addOptional(P,'NewHeader',[],@(x) iscell(x) || ischar(x));
addParameter(P,'PlotOn','ploton',@(x) strcmpi(x,'ploton') || strcmpi(x,'plotoff'));
addParameter(P,'SkipLink','directlink',@(x) strcmpi(x,'directlink') || strcmpi(x,'treelink'));
addParameter(P,'DistMode','shmham',@(x) strcmpi(x,'shmham') || strcmpi(x,'ham'));
addParameter(P,'OutputDistMode','shmham',@(x) strcmpi(x,'shmham') || strcmpi(x,'ham'));
addParameter(P,'SaveName','',@ischar);
addParameter(P,'SavePath','',@ischar);
parse(P,varargin{:});

VDJdata = P.Results.VDJdata;
NewHeader = P.Results.NewHeader;
PlotOn = P.Results.PlotOn;
SkipLink = P.Results.SkipLink;
DistMode = P.Results.DistMode;
OutputDistMode = P.Results.OutputDistMode;
SaveName = P.Results.SaveName;
SavePath = P.Results.SavePath;
if isempty(SaveName)
    [SaveName SavePath] = uiputfile('*.png','Save the tree files as');
end
if isempty(SavePath)
    if ispc
        SavePath = [cd '\'];
    else
        SavePath = [cd '/'];
    end
end

%Select the file or VDJdata
if ischar(VDJdata) || (iscell(VDJdata) && (ischar(NewHeader) || isempty(NewHeader))) %It's not really VDJdata, but a file name.
    FileNames = VDJdata; 
    if ~iscell(FileNames)
        FileNames = {FileNames};
    end
    
    FilePath = NewHeader;
    if isempty(FilePath)
        if ispc
            FilePath = [cd '\'];
        else
            FilePath = [cd '/'];
        end
    end
    SkipLoad = 0;
elseif isempty(VDJdata)
    [FileNames,FilePath] = uigetfile('*.xlsx','multiselect','on');
    if ~iscell(FileNames)
        FileNames = {FileNames};
    end
    SkipLoad = 0;
else
    FileNames = ' ';
    SkipLoad = 1;
end

%==========================================================================
if SkipLoad == 0
    AllVDJdata = cell(1,length(FileNames));
    for f = 1:length(FileNames)    
        [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
        getHeaderVar;

        %Basic data formatting
        VDJdata = removeNAN(VDJdata);
        for j = 1:size(VDJdata,1)
            if isempty(VDJdata{j,TemplateLoc})
                VDJdata{j,TemplateLoc} = 1;
            end
        end
        AllVDJdata{f} = removeDupVDJdata(VDJdata,NewHeader);        
    end
else
    getHeaderVar;
    if iscell(VDJdata{end,end}) %you have a cell in a cell;
        AllVDJdata = VDJdata;
    else
        AllVDJdata = {VDJdata};
    end    
end
clear VDJdata;

%Go through AllVDJdata, to get all unique CDR3s
for f = 1:length(AllVDJdata)
    if f == 1;
        CDR3seq = AllVDJdata{f}(:,CDR3Loc(1));
    else
        CDR3seq = [CDR3seq; AllVDJdata{f}(:,CDR3Loc(1))];
    end
end
[UnqCDR3,~,~] = plotTree_getUnqCDR3(CDR3seq);
[UnqClr,~] = plotTree_getDotClr(UnqCDR3);

%Build the tree, but compute the color scheme based on the unique CDR3 clr
UnqNum = 1; %Unique Picture Number, in case you want to save tree with same CDR3 name.
for f = 1:length(AllVDJdata)
    VDJdata = AllVDJdata{f};    
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);

    if strcmpi(SkipLink,'directlink')
        AncMapCellT = calcAncMapCell(VDJdata,NewHeader,OutputDistMode);

        %Sort by tree height
        SizeMat = zeros(size(AncMapCellT,1),1);
        for k = 1:size(AncMapCellT,1)
             SizeMat(k) = size(AncMapCellT{k},1);
        end
        [SizeMat,SortIdx] = sort(SizeMat,'descend');
        AncMapCellT = AncMapCellT(SortIdx,:);
    end
    
    VDJnames = cell(length(GrpNumUnq),1);
    for y = 1:length(GrpNumUnq)
        IdxLoc = find(GrpNumUnq(y) == GrpNum);    
        Tdata = VDJdata(IdxLoc,:);

        %Include the Group name, VDJ names, and shorten them
        VDJname = Tdata(1,FamLoc);
        for w = 1:3
            TempName = VDJname{w};
            ColonLoc = regexp(TempName,'\:');
            if isempty(ColonLoc)
                ColonLoc = 0;
            end
            TempName = TempName(ColonLoc+1:end);
            TempName = strrep(TempName,' ','');            
            TempName = regexp(TempName,'\|','split');
            TempName = strrep(TempName{1},'IGH','');
            
            %Remove the "0" in V01, D01, etc.
            if TempName(2) == '0'
                TempName(2) = [];
            end
            VDJname{w} = TempName;
        end
        VDJnames{y} = sprintf('%s | %s | %s, Size=%d',VDJname{1},VDJname{2},VDJname{3},size(Tdata,1));
    end
    VDJnames = VDJnames(SortIdx); %need to sort this to match with AncMapCellT ordering
    
    for y = 1:length(GrpNumUnq)
        AncMapCell = AncMapCellT(y,:);
        TreeName = VDJnames{y};
        if strcmpi(PlotOn,'ploton')
            for t = 1:size(AncMapCell,1)
                CDR3seq = AncMapCell{t,2};

                %Assign a DotColor for every entity in AncMap
                DotClr = zeros(size(CDR3seq,1),3);
                for j = 1:size(CDR3seq,1)
                    %Locate the UnqCDR3seq location
                    UnqLoc = 0;
                    for k = 1:size(UnqCDR3,1)
                        if sum(CDR3seq{j}~=UnqCDR3{k}) == 0
                            UnqLoc = k;
                            break
                        end
                    end
                    if UnqLoc > 0
                        DotClr(j,:) = UnqClr(UnqLoc,:);
                    end
                end
                
                [Gx,Ax] = plotTree2(AncMapCell{t,1},TreeName,AncMapCell{t,2},DotClr,'SortMode','sort','FigFontSize',8,'SetXlim',25,'SetYincr',0.01,'FigWidth',2.3,'ScaleFactor',20,'MaxDotArea',400);
                if length(GrpNumUnq) == 0
                    set(Ax,'XTickLabel','');
                    delete(get(Ax,'xlabel'));
                    PlotPos = get(Ax,'Position');
                    PlotPos(2) = 0.014;
                    set(Ax,'Position',PlotPos);
                    
                    TightPos = get(Ax,'TightInset');
                    FigPos = get(Gx,'PaperPosition');
                    FigPos(4) = sum(PlotPos([2 4])) + TightPos(4); %Want to ensure some border
                    FigPos(2) = 0;
                    set(Gx,'Position',FigPos);
                    set(Gx,'PaperPosition',FigPos);
                    set(Gx,'PaperSize',FigPos(3:4));

                elseif strcmpi(OutputDistMode,'shmham') && strcmpi(SkipLink,'directlink')
                    set(get(Ax,'xlabel'),'string','SHM Distance');
                end
                if ~exist('FileName','var')
                    FileName = SaveName;
                    FilePath = SavePath;
                end
                DotLoc = regexp(FileName,'\.');
                saveas(Gx,[FilePath FileName(1:DotLoc(end)-1) '_' upper(DistMode) '_' CDR3seq{1} num2str(UnqNum) '.fig']);
                print(Gx,[FilePath FileName(1:DotLoc(end)-1) '_' upper(DistMode) '_' CDR3seq{1} num2str(UnqNum) '.png'],'-dpng','-r600');
                UnqNum = UnqNum + 1;
                close(Gx);
            end
        end
    end
end