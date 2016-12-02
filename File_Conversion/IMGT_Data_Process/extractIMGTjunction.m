%extractIMGTjunction will look through the
%IMGT_HighV-QUEST_individual_files_folder that has the full Seq output
%data, and then extract the necessary information from the Junction
%Analysis portion that can be used to translate IMGT results to VDJdata
%format. The conversion step is done via the convertIMGT2VDJdata function. 

function extractIMGTjunction()
[Vmap,Dmap,Jmap] = getCurrentDatabase;

%Finding IMGT full sequence files
FileList = dir('*');
DelList = zeros(length(FileList),1)>1;
for j = 1:size(FileList)
    FileName = FileList(j).name;
    if isempty(regexp(FileName,'Seq_','once'))
        DelList(j) = 1;
    end
end
FileList(DelList) = [];

%Perform sorting based on VDJdata sequence number (for ease of matching);
FileOrder = struct2cell(FileList)';
FileMap = zeros(size(FileOrder,1),2); %From this to that
for j = 1:size(FileMap,1)
    SeqOrdT = regexp(FileOrder{j,1},'_','split');
    FileMap(j,1) = j;
    FileMap(j,2) = eval(SeqOrdT{3});
end
FileMap = sortrows(FileMap,2);
FileList = FileList(FileMap(:,1));

%Generate the IMGTdata cell;
IMGTheader = {'SeqName', 'nucleotide', 'V_name', 'D_name', 'J_name', '3''V-REGION', 'Pv', 'N1', 'Pd5', 'D-REGION', 'Pd3', 'N2', 'Pj', '5''J-REGION', 'Vmut', 'Dmut', 'Jmut'};
IMGTdata = cell(length(FileList),length(IMGTheader));
for j = 1:size(FileList)
    FileName = FileList(j).name;
    UndLoc = find(FileName == '_');
    UndLoc = UndLoc(end);
    IMGTdata{j,1} = FileName(1:UndLoc-1); %Samplename is filename
    if ~isempty(regexp(FileName,'Seq','once'))
        FID = fopen(FileName,'r');
        
        %------------------------------------------------------------------
        %Identify Locations of key Data.
            % A. Detailed results  -  get the Sequence Data
            % 1. Alignment for V-GENE  -  get the Vseg + Vref
            % 2. Alignment for D-GENE  -  get the Dseg + Dref + N1/N2
            % 3. Alignment for J-GENE  -  get the Jseg + Jseg
            % 4. Results of IMGT/JunctionAnalysis %Get All
        KeyLocs = zeros(5,1);
        HeaderLabel = {'Sequence number' '1. Align' '2. Align', '3. Align' '4. Results'};
        FoundCt = 0;
        JunctionAnalysis = 1; %1 if JA is done, 0 if JA is not done.
        while feof(FID) == 0
            TextLine= fgetl(FID);            
            for q = 1:5
                if KeyLocs(q) == 0
                    if strfind(TextLine,HeaderLabel{q});
                        KeyLocs(q) = ftell(FID);
                        FoundCt = FoundCt + 1;
                    end
                end
            end
            
            %Check if there is no results on JunctionAnalysis
            if strfind(TextLine,'JunctionAnalysis;No results')
                JunctionAnalysis = 0;
            end
            
            if FoundCt == length(HeaderLabel); break; end %
        end
        if KeyLocs(5) == 0
            JunctionAnalysis = 0; %Sometimes, IMGT will shows JA but then say no JA can be done. Sometimes it will just not even have this section.
        end
        
        %------------------------------------------------------------------        
        %Extract the sample name and sequence in the FASTA format
        fseek(FID,KeyLocs(1),'bof');
        InsFound = 0; %Presence of insertion in sequence, marked as capital in IMGT Fasta Seq format
        DelFound = 0;
        while ftell(FID) <= KeyLocs(2)
            TextLine = fgetl(FID);
            % Extract sequence name and sequence from FASTA fromat
            if ~isempty(strfind(TextLine,'>'))
                SampName = TextLine(2:end);
                SampSeq = '';
                while ~isempty(TextLine) %IMGT leaves blanks to separate paragraphs
                    TextLine= fgetl(FID);
                    SampSeq = [SampSeq TextLine];
                end
                
                %Look for insertions here, marked as capital.
                InsLoc = isstrprop(SampSeq,'upper');
                InsLen = sum(InsLoc);
                if sum(InsLen) > 0
                    InsFound = 1; %triggers the full V alignment search;
                end
                
                %Save variables to IMGTdata
                SampSeq = upper(SampSeq);
%                IMGTdata{j,1} = SampName; Changed this so that file name
%                is seqname
                IMGTdata{j,2} = SampSeq;
                
                break                
            end
        end
        while ftell(FID) <= KeyLocs(2)
            TextLine = fgetl(FID);            
            %Look for deletion information in the results summary
            if ~isempty(strfind(TextLine,'deletions have'))
                DelFound = 1;
                break
            end
        end
        
        %------------------------------------------------------------------        
        %If InDel is found, go through the Valignement to correct InDel
        if DelFound == 1 || InsFound == 1
            fseek(FID,KeyLocs(2),'bof');
            
            %Extract the alignment information
            CutLoc = 0;
            NewSeq = '';
            RefSeq = '';
            while feof(FID) == 0
                TextLine= fgetl(FID);
                
                %Set the hard while stops
                if sum(regexp(TextLine,'\d+\.\s+','once')==1) > 0; break; end
                                
                if ~isempty(strfind(TextLine,SampName))
                    %Determine IMGT format cut location (do this once)
                    if CutLoc == 0;
                        TextLine = strrep(TextLine,SampName,repmat(' ',1,length(SampName)));
                        CutLoc = regexp(TextLine,'[^\s]','once');
                    end
                    
                    %Extract the sequence and reference alignment info
                    SeqEval = TextLine(CutLoc:end);
                    NewSeq = [NewSeq SeqEval];
                    TextLine= fgetl(FID); %The best ref seq is just next one.
                    RefSeq = [RefSeq TextLine(CutLoc:end)];
                end
            end
            
            %Trim the NewSeq, replace the . with x, and trim the length to conserve original length.
            StartLoc = regexp(NewSeq,'\w','once') - InsLen; %Will have to "steal" nts from the ref seq
            NewSeq = NewSeq(StartLoc:end);
            RefSeq = RefSeq(StartLoc:end);
            
            %Make sure to substitute the inserts before modifying '.'
            NewSeq(1:InsLen) = RefSeq(1:InsLen);
            
            %Remove the "..." dots from the IMGT standard notation
            RefDotLoc = regexp(RefSeq,'\.');
            NewSeq(RefDotLoc) = [];
            RefSeq(RefDotLoc) = [];
            
            %Replace the seqeunce deletions with the expect NT from RefS
            NewDotLoc = regexp(NewSeq,'\.');
            NewSeq(NewDotLoc) = RefSeq(NewDotLoc);

            IMGTdata{j,2} = upper(NewSeq);
        end
            
        %------------------------------------------------------------------
        %If no Junction Analysis, get V,D,J alignment info
        if JunctionAnalysis == 0
            if KeyLocs(2) > 0
                [Vseq, Vref, Vlen, Vdel5, Vdel3, Vname] = getIMGTalign(FID,SampName,KeyLocs(2));
            else
                Vseq = '';
                Vref = '';
                Vlen = 0;
                Vdel3 = 0;
                Vdel5 = 0;
                Vname = 'Unresolved';
            end
            if KeyLocs(3) > 0
                [Dseq, Dref, Dlen, Ddel5, Ddel3, Dname] = getIMGTalign(FID,SampName,KeyLocs(3));
                Dref(Dref=='-') = Dseq(Dref=='-');
                
                [DrefSeq,DmapNum] = lookupVDJgene(Dname,Vmap,Dmap,Jmap);
                %Determine the true Ddel3,Ddel5
                Dmatch = findGeneMatch(Dref,Dmap(DmapNum,:),'D');
                Ddel5 = Dmatch{3}(1);
                Ddel3 = Dmatch{3}(3);
                Dlen = Dmatch{3}(2);
            else
                Dseq = '';
                Dref = '';
                Dlen = 0;
                Ddel3 = 0;
                Ddel5 = 0;
                Dname = 'Unresolved';
            end
            if KeyLocs(4) > 0
                [Jseq, Jref, Jlen, Jdel5, Jdel3, Jname] = getIMGTalign(FID,SampName,KeyLocs(4));                        
            else
                Jseq = '';
                Jref = '';
                Jlen = 0;
                Jdel3 = 0;
                Jdel5 = 0;
                Jname = 'Unresolved';
            end
        
            VMDNJ = [Vlen 0 0 0 Jlen];
            
            if sum(VMDNJ) > length(SampSeq);
                %Maybe V and J have too many nts mismatched near the ends
                DashLoc = find(Vref == '-');
                DiffLoc = find(diff(DashLoc)>= 4);
                if ~isempty(DiffLoc)
                    EndLoc = DashLoc(DiffLoc(end))-1;
                    Vseq = Vseq(1:EndLoc);
                    Vref = Vref(1:EndLoc);
                    Vlen = length(Vseq);
                else
                    continue %Something is wrong
                end
                VMDNJ = [Vlen 0 0 0 Jlen];
            end
            
            while (length(SampSeq) - sum(VMDNJ)) <= length(Dseq)
                MissLoc = regexpi(Vref,'[^\-]');
                if isempty(MissLoc); break; end
                for q = length(MissLoc):-1:2
                    if MissLoc(q) == MissLoc(q-1)+1
                        MissLoc(q) = [];
                    else
                        MissLoc = MissLoc(q);
                        break
                    end
                end
                AddDel = length(Vref) - MissLoc + 1;
                Vdel3 = Vdel3 + AddDel;
                VMDNJ(1) = VMDNJ(1) - AddDel;
            end
               
            if (sum(VMDNJ([1 end])) <= length(SampSeq)) && (KeyLocs(2) * KeyLocs(3) * KeyLocs(4)) > 0
                MDNnt = SampSeq(VMDNJ(1)+1:end-VMDNJ(end));
                [~,~,~,MatchAt] = convolveSeq(MDNnt,Dseq,0,0);
                Mlen = MatchAt(1) - 1;
                Nlen = length(MDNnt) - MatchAt(2);
                VMDNJ(2:4) = [Mlen length(Dseq) Nlen];
            end
                        
            %Fill in IMGT table
            IMGTdata(j,3:5) = {Vname Dname Jname};      
            
            if min(VMDNJ) < 0 || sum(VMDNJ) ~= length(SampSeq)
                fclose(FID);
                continue %Something is wrong
            end
            
            IMGTdata{j,6} = repmat('.',1,Vdel3); %Vinfo (this does not include all nts.
            
            if VMDNJ(2) > 0 %N2 region
                IMGTdata{j,8} = MDNnt(1:VMDNJ(2));
            else
                IMGTdata{j,8} = '';
            end
            
            IMGTdata{j,10} = [repmat('.',1,Ddel5) SampSeq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))) repmat('.',1,Ddel3)]; %D region
            
            if VMDNJ(4) > 0 %N1 region
                IMGTdata{j,12} = MDNnt(end-VMDNJ(4)+1:end); 
            end
            
            IMGTdata{j,14} = [repmat('.',1,Jdel5) SampSeq(end-VMDNJ(end)+1:end)]; %
        end
        
        if JunctionAnalysis == 1
            fseek(FID,KeyLocs(5),'bof');

            %Extract the VDJ genes from Junction Analysis
            TopLine = '';
            BotLine = '';
            CutLoc = 0;
            while feof(FID) == 0
                TextLine= fgetl(FID);

                %Set the while stops
                if sum(regexp(TextLine,'\d+\.\s+','once')==1) > 0; break; end

                %Collect the IMGT formatted lines into 2 lines
                if ~isempty(strfind(TextLine,'Input'))
                    if CutLoc == 0;
                        TextLine = strrep(TextLine,'Input','     ');
                        CutLoc = regexp(TextLine,'[^\s]','once');                        
                    end

                    TopLine = [TopLine TextLine(CutLoc:end)];
                    TextLine= fgetl(FID);                       
                    BotLine = [BotLine TextLine(CutLoc:end)];
                end
            end

            %Replace single space fields with _ to prevent bad splits
            SingSpaceLoc = regexp(TopLine,'[\w\d]\s{1,1}[\w\d]+');
            TopLine(SingSpaceLoc+1) = '_';                        
            SingSpaceLoc = regexp(BotLine,'[\w\d]\s{1,1}[\w\d]+');
            BotLine(SingSpaceLoc+1) = '_';

            %Split the text to extract the data
            FieldName = regexp(TopLine,'\s+','split');
            FieldData = regexp(BotLine,'\s*','split');

            %Fill the the data table according to field names
            for k = 1:length(FieldName)
                %Make unique P labels, which is not done by default
                if strcmpi(FieldName{k},'P')
                    if strcmpi(FieldName{k-1},'3''V-REGION')
                        FieldName{k} = 'Pv';
                    elseif strcmpi(FieldName{k-1},'N1')
                        FieldName{k} = 'Pd5';
                    elseif strcmpi(FieldName{k-1},'D-REGION')
                        FieldName{k} = 'Pd3';
                    elseif strcmpi(FieldName{k-1},'N2')
                        FieldName{k} = 'Pj';
                    end
                end                    

                FieldLoc = findHeader(IMGTheader,FieldName{k});
                if FieldLoc == 0; continue; end %Non existent field
                IMGTdata{j,FieldLoc} = FieldData{k};
            end
        end        
        fclose(FID);
        j
    end
end


%Save the extracted data as a file.
SaveName = 'ExtractedJunctionData';
if ispc
    xlswrite([SaveName '.xlsx'],cat(1,IMGTheader,IMGTdata));
else
    writeDlmFile(cat(1,IMGTheader,IMGTdata),[SaveName '.csv'],'\t');
end