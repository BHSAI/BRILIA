%buildTree_SameColorCode will go through VDJdata and plot the trees per unique group.
%HOWEVER, but will first look for all unique sequences in both files,
%assemble the color code scheme, and then apply them throughout the entire
%trees for all files. This is useful if you have a main tree, subtrees from
%2 files, and you want to color code the CDR3 regions the same way.

%1) Identify lowest mutated V and J and D.
%2) Assemble reference sequence
%3) Calculated distances from ref to seq
%4) Calculated distances from seq to seq1
%etc etc

function VDJdata = buildTree_SameColorCode(varargin)
%Check input for options
PlotOn = 'ploton';
SaveOn = 'saveoff';
DistMode = 'shmham';
SkipLink = 'directlink'; %skiplink will use buildTreeLink, directlink will use the direct ref-seq pairing of VDJdata.
OutputDistMode = 'shmham';
DelThis = zeros(1,length(varargin)) > 1;
if ~isempty(varargin)
    for k = 1:length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'ploton'
                    PlotOn = 'ploton';
                    DelThis(k) = 1;
                case 'plotoff'
                    PlotOn = 'plotoff';
                    DelThis(k) = 1;
                case 'saveon'
                    SaveOn = 'saveon';
                    DelThis(k) = 1
                case 'saveoff'
                    SaveOn = 'saveoff';
                    DelThis(k) = 1;
                case 'ham'
                    DistMode = 'ham';
                    DelThis(k) = 1;
                case 'shm'
                    DistMode = 'shm';
                    DelThis(k) = 1;
                case 'shmham'
                    DistMode = 'shmham';
                    DelThis(k) = 1;
                case 'none'
                    DistMode = 'none';
                    DelThis(k) = 1;
                case 'skiplink'
                    SkipLink = 'treelink';
                    DelThis(k) = 1;
                case 'directlink'
                    SkipLink = 'directlink';
                    DelThis(k) = 1;
            end
        end
    end
end
Inputs = varargin(DelThis == 0); %Leftover inputs

%Select the file or VDJdata
if length(Inputs)==2
    VDJdata = Inputs{1};
    NewHeader = Inputs{2};
    SkipLoad = 1;
    FileNames = {' '};
else
    if length(varargin) == 1 %You've listed the file anmes'
        FileNames = varargin{1}{1};
        FilePath = '';
        SkipLoad = -1;
    else %open the files names
        [FileNames,FilePath] = uigetfile('*.xlsx','multiselect','on');
        SkipLoad = 0;
    end
    if ~iscell(FileNames)
        FileNames = {FileNames};
    end
end

%Pick the bigger tree file
for k = 1:length(FileNames)
    disp([num2str(k) ') ' FileNames{k}]);
end
FirstFile = input('Which one is the bigger tree file? ');
if FirstFile == 2
    FileNames = FileNames([FirstFile 1]);
end

%==========================================================================
for f = 1:length(FileNames)    
    if SkipLoad <= 0
        if SkipLoad < 0 %It's included in files names
            [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FileNames{f}]);
        else
            [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
        end
    end
    getHeaderVar;
    
    %Basic data formatting
    VDJdata = removeNAN(VDJdata);
    for j = 1:size(VDJdata,1)
        if isempty(VDJdata{j,TemplateLoc})
            VDJdata{j,TemplateLoc} = 1;
        end
    end
    VDJdata = removeDupVDJdata(VDJdata,NewHeader);
    
    if f == 1;
        CDR3seq = VDJdata(:,CDR3Loc(1));
    else
        CDR3seq = [CDR3seq; VDJdata(:,CDR3Loc(1))];
    end
end

[UnqCDR3,~,~] = plotTree_getUnqCDR3(CDR3seq);
[UnqClr,~] = plotTree_getDotClr(UnqCDR3);

%Build the tree, but compute the color scheme based on the unique CDR3 clr
for f = 1:length(FileNames)    
    if SkipLoad == 0
        [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
    end
    getHeaderVar;
    
    %Basic data formatting
    VDJdata = removeNAN(VDJdata);
    for j = 1:size(VDJdata,1)
        if isempty(VDJdata{j,TemplateLoc})
            VDJdata{j,TemplateLoc} = 1;
        end
    end
    VDJdata = removeDupVDJdata(VDJdata,NewHeader);

    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);
    UnqNum = 1; %Unique Picture Number, in case you want to save tree with same CDR3 name.
    if strcmpi(SkipLink,'directlink')
        AncMapCellT = calcAncMapCell(VDJdata,NewHeader,OutputDistMode); %USING SHM distance!

        %Sort by tree height, and then, figure out which one will have the
        %X labels.
        SizeMat = zeros(size(AncMapCellT,1),1);
        for k = 1:size(AncMapCellT,1)
             SizeMat(k) = size(AncMapCellT{k},1);
        end
        [SizeMat,SortIdx] = sort(SizeMat,'descend');
        AncMapCellT = AncMapCellT(SortIdx,:);

        DivideLoc = 3; %Determined by hand.
        if size(AncMapCellT,1) == 1
            DivideLoc = 1;
        end
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
            DotLoc = regexp(FileName,'\.');
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
                if (y == DivideLoc || y == length(GrpNumUnq)) == 0
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
                saveas(Gx,[FilePath FileName(1:DotLoc(end)-1) '_' upper(DistMode) '_' CDR3seq{1} num2str(UnqNum) '.fig']);
%                print(Gx,[FilePath FileName(1:DotLoc(end)-1) '_' upper(DistMode) '_' CDR3seq{1} num2str(UnqNum) '.png'],'-dpng','-r600');
                UnqNum = UnqNum + 1;
                close(Gx);
            end
        end
    end
end