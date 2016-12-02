%buildTree will go through VDJdata,

%1) Identify lowest mutated V and J and D.
%2) Assemble reference sequence
%3) Calculated distances from ref to seq
%4) Calculated distances from seq to seq1
%etc etc

function VDJdata = buildTree(varargin)
%Check input for options
PlotOn = 'ploton';
SaveOn = 'saveoff';
DistMode = 'shmham';
SkipLink = 'directlink'; %directlink will use the VDJdaata RefSeq-Seq relation to build the trees.
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
                    DelThis(k) = 1;
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
    [FileNames,FilePath] = uigetfile('*.xlsx','multiselect','on');
    if ~iscell(FileNames)
        FileNames = {FileNames};
    end
    SkipLoad = 0;
end

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
    
%     %Remove nonfuncitonal VDJs
%     NonprodLoc = char(VDJdata(:,FunctLoc)) == 'N';
%     VDJdata = VDJdata(NonprodLoc==0,:); %productives only
    
    %Optional VDJdata corrections
    RenumberThis = 0;
    if RenumberThis == 1
        VDJdata = renumberVDJmap(VDJdata,NewHeader); %Just in case the database changed. Need to renumber the reference map number.
        VDJdata = rankVDJdata(VDJdata,NewHeader); %This should already be ranked due to the conformGeneGroup function.
        VDJdata(:,RefSeqLoc) = buildRefSeq(VDJdata,NewHeader);
    end

    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);
    UnqNum = 1; %Unique Picture Number, in case you want to save tree with same CDR3 name.

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
        TreeName = sprintf('%s | %s | %s, Size=%d',VDJname{1},VDJname{2},VDJname{3},size(Tdata,1));
        
        if strcmpi(SkipLink,'treelink')
            DevPerc = 0.03;
            MaxDev = ceil(length(Tdata{1,SeqLoc})*DevPerc); 
            [AncMapCell,Tdata] = buildTreeLink(Tdata,NewHeader,MaxDev);

            %Go ahead and switch off the ref seq.
            for v = 1:size(AncMapCell,1)
                ChildParMap = AncMapCell{v,1};
                %Update original Tdata based on this ChildPar relationship
                for g = 1:size(ChildParMap,1)
                    if ChildParMap(g,2) > 0
                        Tdata(ChildParMap(g,1),RefSeqLoc) = Tdata(ChildParMap(g,2),SeqLoc);
                    end
                end
            end
            VDJdata(IdxLoc,:) = Tdata;
        else strcmpi(SkipLink,'directlink')
            AncMapCell = calcAncMapCell(Tdata,NewHeader);
            for k = 1:size(AncMapCell,1)
                NewAncMap = AncMapCell{k};
                TreeMapTemp= findTreeClust(NewAncMap);
                if max(TreeMapTemp(:,2)) > 1; 
                    error('You have multiple clusters')
                elseif sum(NewAncMap(:,1) == NewAncMap(:,2)) > 0
                    error('parent child same')
                elseif sum(NewAncMap(:,2) == 0 ) > 1
                    error('multiple roots')
                end
            end
        end
           
        if strcmpi(PlotOn,'ploton')
            DotLoc = regexp(FileName,'\.');
            for t = 1:size(AncMapCell,1)
                CDR3seq = AncMapCell{t,2};
                if size(CDR3seq,1) < 2; continue; end
                [Gx,Ax] = plotTree2(AncMapCell{t,1},TreeName,AncMapCell{t,2},'SortMode','sort');
                print(Gx,[FilePath FileName(1:DotLoc(end)-1) '_' upper(DistMode) '_' CDR3seq{1} num2str(UnqNum) '.png'],'-dpng','-r600');
                UnqNum = UnqNum + 1;
                close(Gx);
            end
%             %Remember that if you use directlink, you only have to do this
%             %once. End the for loop to do all groups.
%             if strcmpi(SkipLink,'directlink')
%                 break
%             end
        end
    end
    
    if strcmpi(SaveOn,'saveon')
        %Before saving to xlsx, convert columns with matrix values into char
        VDJdata = reformatAlignment(VDJdata,1);
        for q = 1:size(VDJdata,1)
            for w = 1:3
                VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
            end
        end
        if SkipLoad == 0
            DotLoc = regexp(FileName,'\.');
            xlswrite([FilePath FileName(1:DotLoc(end)-1) 'tree.xlsx'],[NewHeader; VDJdata]);
        else
            [FileName,FilePath] = uigetfile('*.xlsx','Save the tree excel sheet as');
            DotLoc = regexp(FileName,'\.');
            xlswrite([FilePath FileName(1:DotLoc(end)-1) 'tree.xlsx'],[NewHeader; VDJdata]);
        end
    end
end