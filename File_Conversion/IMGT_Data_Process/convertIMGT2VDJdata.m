%convertIMGT2VDJdata will take the ExtractedJunctionData.xlsx/csv file and
%then convert it VDJdata format. Before running this, run the
%extractJunctionData function on the IMGT HighV-QUEST results. 

function convertIMGT2VDJdata()
%Open the Junction Data
[FileName, FilePath] = uigetfile('*.xlsx;*.csv','Open the IMGT junction file','ExtractedJunctionData');
[IMGTdata, IMGTheader, ~, ~] = openSeqData([FilePath FileName]);
for j = 1:size(IMGTdata,1)
    for k = 1:size(IMGTdata,2)
        if isnan(IMGTdata{j,k})
            IMGTdata{j,k} = '';
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

%IMGTdata column locations
SeqNameLocI = findHeader(IMGTheader,'SeqName');
SeqLocI = findHeader(IMGTheader,'nucleotide');
FamLocI = findHeader(IMGTheader,{'V_name','D_name','J_name'});
VsegLoc = findHeader(IMGTheader,'3''V-REGION');
MsegLoc = findHeader(IMGTheader,{'Pv' 'N1' 'Pd5'});
DsegLoc = findHeader(IMGTheader,'D-REGION');
NsegLoc = findHeader(IMGTheader,{'Pd3' 'N2' 'Pj'});
JsegLoc = findHeader(IMGTheader,'5''J-REGION');

%Build up the VDJdata matrix
VDJdata = cell(size(IMGTdata,1),length(NewHeader));
KeepThese = ones(size(IMGTdata,1),1) == 1;
for j = 1:size(IMGTdata,1)
    FullSeq = IMGTdata{j,SeqLocI};
    VDJdata{j,SeqNumLoc} = eval(cell2mat(regexp(IMGTdata{j,SeqNameLocI},'[\d*]+','match')));
    VDJdata(j,SeqLoc) = IMGTdata(j,SeqLocI);

    %Determine gene names and number
    Vname = 'Unresolved';
    VmapNum = 0;
    if ~isempty(IMGTdata{j,FamLocI(1)})        
        NameSeg = regexp(IMGTdata{j,FamLocI(1)},'\_|\s','split');
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
    if ~isempty(IMGTdata{j,FamLocI(2)})
        NameSeg = regexp(IMGTdata{j,FamLocI(2)},'\_|\s','split');
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
    if ~isempty(IMGTdata{j,FamLocI(3)})
        NameSeg = regexp(IMGTdata{j,FamLocI(3)},'\_|\s','split');
        for q = 1:length(NameSeg)
            if sum(regexpi(NameSeg{q},'IGHJ')) > 0
                Jname = NameSeg{q};
                break
            end
        end
        
        if ~strcmpi(Jname,'unresolved')
            [~, JmapNum] = lookupVDJgene(Jname,Vmap,Dmap,Jmap);
            Jname = Jmap{JmapNum,3};
        end
    end
    VDJdata{j,FamNumLoc(3)} = JmapNum;
    VDJdata{j,FamLoc(3)} = Jname;

    if VmapNum(1) * JmapNum(1) * DmapNum(1) == 0  %Unresolved or unknown family cases, skip.
        KeepThese(j) = 0;
        continue
    end

    %Determine the VMDNJ segment lengths
    Vnts = cell2mat(IMGTdata(j,VsegLoc));
    Mnts = cell2mat(IMGTdata(j,MsegLoc));
    Dnts = cell2mat(IMGTdata(j,DsegLoc)); 
    Nnts = cell2mat(IMGTdata(j,NsegLoc)); 
    Jnts = cell2mat(IMGTdata(j,JsegLoc));
    
    if isempty(Vnts) || isempty(Dnts) || isempty(Jnts) %No junction found.
        KeepThese(j) = 0;
        continue;
    end
    
    %Determine the deletion counts, marked as "." 
    Vdel = 0;
    for k = length(Vnts):-1:1
        if Vnts(k) == '.'
            Vdel = Vdel+1;
        else
            break
        end
    end
    %Always confirm if Vdel is correct, as JA does not have right database.
    VCheck = Vmap{VmapNum}(end-Vmap{VmapNum,end}+1:end);
    FromJA = 1;
    if sum(Vnts ~='.') == 0 %Then this is NOT from JA, but form VQUEST
        FromJA = 0;
    end
    if length(VCheck) ~= length(Vnts) && FromJA == 1
        Vdel = Vdel + length(VCheck) - length(Vnts);
    end
    
    
    Ddel5 = 0;
    for k = 1:length(Dnts)
        if Dnts(k) == '.'
            Ddel5 = Ddel5+1;
        else
            break
        end
    end
    
    Ddel3 = 0;
    for k = length(Dnts):-1:1
        if Dnts(k) == '.'
            Ddel3 = Ddel3+1;
        else
            break
        end
    end
    
    %Make sure Del5+D+Del3 == DalignLen
    DrefSeq = Dmap{DmapNum,1};
    DrefLen = length(DrefSeq);
    DsamLen = Ddel3+Ddel5+sum(Dnts ~= '.');
    if DrefLen ~= DsamLen %Recalc the deletion edges
        Dmatch = findGeneMatch(Dnts(Dnts ~= '.'),Dmap(DmapNum,:),'D');
        Ddel5 = Dmatch{3}(1);
        Ddel3 = Dmatch{3}(3);
        Dnts = [repmat('.',1,Ddel5) Dnts(Dnts ~= '.') repmat('.',1,Ddel3)];
    end
    
    Jdel = 0;
    for k = 1:length(Jnts)
        if Jnts(k) == '.'
            Jdel = Jdel+1;
        else
            break
        end
    end
    %Always confirm if Vdel is correct, as JA does not have right database.
    JCheck = Jmap{JmapNum}(1:Jmap{JmapNum,end}+2);
    if length(JCheck) ~= length(Jnts) && FromJA == 1
        Jdel = Jdel + length(JCheck) - length(Jnts);
    end    
    DelCt = [Vdel Ddel5 Ddel3 Jdel];
    VDJdata(j,DelLoc) = num2cell(DelCt);
    
    Vnts = FullSeq(1:length(FullSeq)-length(Jnts)-length(Nnts)-length(Dnts)-length(Mnts)+Vdel+Jdel+Ddel5+Ddel3);
    VMDNJ = [length(Vnts)-Vdel length(Mnts) length(Dnts)-Ddel5-Ddel3 length(Nnts) length(Jnts)-Jdel];
    VDJdata(j,LengthLoc) = num2cell(VMDNJ);
    
    if isempty(Jnts) || isempty(Dnts) || isempty(Vnts)
        KeepThese(j) = 0;
        continue
    end
    
end

%Set group number same as seq number
if ~isnumeric(VDJdata{1,SeqNumLoc});
    for j = 1:size(VDJdata,1)
        if isempty(VDJdata{j,SeqNumLoc});
            continue
        end
        VDJdata{j,SeqNumLoc} = eval(VDJdata{j,SeqNumLoc});
    end
end
VDJdata(:,GrpNumLoc) = VDJdata(:,SeqNumLoc);
VDJdata(KeepThese,:) = updateVDJdata(VDJdata(KeepThese,:),NewHeader,Vmap,Dmap,Jmap);
VDJdata = sortrows(VDJdata,SeqNumLoc);

%Save the files
DotLoc = find(FileName == '.');
DotLoc = DotLoc(end);
SaveName = FileName(1:DotLoc-1);

%Before saving to xlsx, convert columns with matrix values into char
save([FilePath SaveName '.IMGT.mat'],'VDJdata','NewHeader')
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
    end
end
if ispc
    xlswrite([FilePath SaveName '.IMGT.xlsx'],[NewHeader; VDJdata]);
else
    writeDlmFile([NewHeader;VDJdata],[FilePath SaveName '.IMGT.csv'],'\t');
end