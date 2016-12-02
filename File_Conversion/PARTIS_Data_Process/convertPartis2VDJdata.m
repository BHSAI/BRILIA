%convertPartis2VDJdata will convert Partis output format to VDJdata format.
%This requires the HeadersMap.xlsx file to find relationship between the
%datasets.

%convertIMGT2VDJdata will take the ExtractedJunctionData.xlsx/csv file and
%then convert it VDJdata format. Before running this, run the
%extractJunctionData function on the IMGT HighV-QUEST results. 

function convertPartis2VDJdata()
%Open the partis Data File (converted to xlsx file for ease).
[FileName, FilePath] = uigetfile('*.xlsx;*.csv','Open the IMGT junction file','ExtractedJunctionData');
[PartisData, PartisHeader, ~, ~] = openSeqData([FilePath FileName]);
for j = 1:size(PartisData,1)
    for k = 1:size(PartisData,2)
        if isnan(PartisData{j,k})
            PartisData{j,k} = '';
        end
    end
end

%Open the VDJ Gene database in .mat file, if it exist. 
[Vmap, Dmap, Jmap] = getCurrentDatabase('change'); %Selective active full database

%Extract the VDJdata headers
[~, ~, StandardData] = xlsread('Headers_Brilia.xlsx');
NewHeaderLoc = findHeader(StandardData(1,:),'VDJdata');
NewHeader = StandardData(2:end,NewHeaderLoc)';
for j = 1:length(NewHeader)
    if isnan(NewHeader{j}); break; end
end
NewHeader(j:end) = [];
getHeaderVar;

%PartisData column locations
SeqNameLocI = findHeader(PartisHeader,'unique_ids');
SeqLocI = findHeader(PartisHeader,'seqs');
RefSeqLocI = findHeader(PartisHeader,'naive_seq');
FamLocI = findHeader(PartisHeader,{'v_gene','d_gene','j_gene'});
VdelLocI = findHeader(PartisHeader,{'v_5p_del','v_3p_del'});
DdelLocI = findHeader(PartisHeader,{'d_5p_del','d_3p_del'});
JdelLocI = findHeader(PartisHeader,{'j_5p_del','j_3p_del'});
NvdLocI = findHeader(PartisHeader,'vd_insertion');
NdjLocI = findHeader(PartisHeader,'dj_insertion');

%Build up the VDJdata matrix
VDJdata = cell(size(PartisData,1),length(NewHeader));
KeepThese = ones(size(PartisData,1),1,'logical');
for j = 1:size(PartisData,1)
    %Fill in the basic seq name and number info
    VDJdata{j,SeqNumLoc} = cell2mat(regexp(PartisData{j,SeqNameLocI},'[\d*]+','match'));
    VDJdata(j,SeqLoc) = PartisData(j,SeqLocI);
    VDJdata(j,SeqNameLoc) = PartisData(j,SeqNameLocI);
    VDJdata(j,RefSeqLoc) = PartisData(j,RefSeqLocI);
    
    %Fill in the VDJ gene name and number info
    Vname = 'Unresolved';
    VmapNum = 0;
    if ~isempty(PartisData{j,FamLocI(1)})        
        NameSeg = regexp(PartisData{j,FamLocI(1)},'\_|\s','split');
        for q = 1:length(NameSeg)
            if sum(regexpi(NameSeg{q},'IGHV')) > 0
                Vname = NameSeg{q};
                break
            end
        end
        
        if ~strcmpi(Vname,'unresolved')
            [~, VmapNum] = lookupVDJgene(Vname,Vmap,Dmap,Jmap);
            Vname = Vmap{VmapNum,3};
        end
    end
    VDJdata{j,FamNumLoc(1)} = VmapNum;
    VDJdata{j,FamLoc(1)} = Vname;

    Dname = 'Unresolved';
    DmapNum = 0;
    if ~isempty(PartisData{j,FamLocI(2)})
        NameSeg = regexp(PartisData{j,FamLocI(2)},'\_|\s','split');
        for q = 1:length(NameSeg)
            if sum(regexpi(NameSeg{q},'IGHD')) > 0
                Dname = NameSeg{q};
                break
            end
        end
        
        if ~strcmpi(Dname,'unresolved')
            [~, DmapNum] = lookupVDJgene(Dname,Vmap,Dmap,Jmap);
            Dname = Dmap{DmapNum,3};
        end
    end
    VDJdata{j,FamNumLoc(2)} = DmapNum;
    VDJdata{j,FamLoc(2)} = Dname;
    
    Jname = 'Unresolved';
    JmapNum = 0;
    if ~isempty(PartisData{j,FamLocI(3)})
        NameSeg = regexp(PartisData{j,FamLocI(3)},'\_|\s','split');
        for q = 1:length(NameSeg)
            if sum(regexpi(NameSeg{q},'IGHJ')) > 0
                Jname = NameSeg{q};
                break
            end
        end
        
        if ~strcmpi(Jname,'unresolved')
            Jname = strrep(Jname,'P',''); %For some reason, partis has a "P" at the J names.
            [~, JmapNum] = lookupVDJgene(Jname,Vmap,Dmap,Jmap);
            Jname = Jmap{JmapNum,3};
        end
    end
    VDJdata{j,FamNumLoc(3)} = JmapNum;
    VDJdata{j,FamLoc(3)} = Jname;
    
    if VmapNum * JmapNum * DmapNum == 0
        KeepThese(j) = 0;
        continue
    end

    %Fill in the VMDNJ length info, and Del info
    Mseq = PartisData{j,NvdLocI};
    Mlen = length(Mseq);
    
    Nseq = PartisData{j,NdjLocI};
    Nlen = length(Nseq);
    
    V5del = PartisData{j,VdelLocI(1)};
    V3del = PartisData{j,VdelLocI(2)};
    VrefLen = length(Vmap{VmapNum,1});
    Vlen = VrefLen - V5del - V3del;
    
    D5del = PartisData{j,DdelLocI(1)};
    D3del = PartisData{j,DdelLocI(2)};
    DrefLen = length(Dmap{DmapNum,1});
    Dlen = DrefLen - D5del - D3del;
    
    J5del =  PartisData{j,JdelLocI(1)};
    J3del =  PartisData{j,JdelLocI(2)};
    JrefLen = length(Jmap{JmapNum,1});
    Jlen = JrefLen - J5del - J3del;
    
    VDJdata(j,LengthLoc) = num2cell([Vlen Mlen Dlen Nlen Jlen]);
    VDJdata(j,DelLoc) = num2cell([V3del D5del D3del J5del]);

end

%Set group number same as seq number
if ~isnumeric(VDJdata{1,SeqNumLoc});
    for j = 1:size(VDJdata,1)
        VDJdata{j,SeqNumLoc} = eval(VDJdata{j,SeqNumLoc});
    end
end
VDJdata(:,GrpNumLoc) = VDJdata(:,SeqNumLoc);

VDJdataT = VDJdata(KeepThese==0,:);
VDJdata = VDJdata(KeepThese,:);
VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 
VDJdata = labelNonprodVDJ(VDJdata,NewHeader);
VDJdata = [VDJdata; VDJdataT]; %Rejoin the neglected ones

VDJdata = sortrows(VDJdata,SeqNumLoc); %Sort it

%Save the files
DotLoc = find(FileName == '.');
DotLoc = DotLoc(end);
SaveName = FileName(1:DotLoc-1);

%Before saving to xlsx, convert columns with matrix values into char
save([FilePath SaveName '.partis.mat'],'VDJdata','NewHeader')
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
    end
end
if ispc
    xlswrite([FilePath SaveName '.partis.xlsx'],[NewHeader; VDJdata]);
else
    writeDlmFile([NewHeader;VDJdata],[FilePath SaveName '.partis.csv'],'\t');
end