%generateSHMlib will take a repertoire, and then introduce SHM according to
%ADAR and AID-mediate mutation patterns. 
%
%  [VDJdata, NewHeader] = generateVDJrep(GrpSize,SHMrate)
%     RepSize = number of unique clones
%     GrpSize = number of Seq per clone
%     SHMrate = number of SHM mutation per 100 bp
%     TDTon = enable/disable TDT insertions. 1 = enable, 0 = disable.
%function [VDJdata, NewHeader] = generateVDJlib(RepSize,GrpSize,SHMrate,TDTon)
function [VDJdata, NewHeader] = generateSHMlib
GrpSize = 5; %Number of sequence per group
SHMrate = 5; %Number of nts to mutate per round
AllowSameMutPos = 0; %0 means do NOT allow same nt to mutate more than 1, or have non-repeating mutaiton locations.
[VDJdata, NewHeader, FileName, FilePath] = openSeqData;

getHeaderVar;

%Probablity nt goes from (column) nt to (row) nt. Based of A Collins 2015
%data set, analyzed by MAVRIC Dev 9.
% Pnt2nt = [...
% 0       10395   25057   9981
% 12809   0       9393    24737
% 34286   10005   0       7375
% 15368   30254   5561    0   ];
% Pnt2nt = Pnt2nt./sum(Pnt2nt(:)); %MAVRIC v9

Pnt2nt = [...
0       7320   22093   5798
9228    0      8200    19160
30650   7874   0       7492
12086   23588  5886    0    ];
Pnt2nt = Pnt2nt./sum(Pnt2nt(:)); %BRILIA v13

Pnt = sum(Pnt2nt,1); %Probablity of nt compositoin mutation
CPnt2nt = cumsum(Pnt2nt,1)./repmat(sum(Pnt2nt,1),4,1); %Probability for pairwise mutation

%--------------------------------------------------------------------------
NewVDJdata = cell(size(VDJdata,1)*GrpSize,size(VDJdata,2));
q = 1; %Counter for NewVDJdata
for j = 1:size(VDJdata,1)
    %Extract basics
    Seq = upper(VDJdata{j,SeqLoc});
    SeqNum = VDJdata{j,SeqNumLoc};
    CDR3start = VDJdata{j,CDR3Loc(3)};
    CDR3end = VDJdata{j,CDR3Loc(4)};
    ReadFrame = mod(CDR3end,3)+1;
    
    MutLoc = zeros(size(Seq)); %Keeps track of where mutations were done
    for g = 1:GrpSize
        %If you allow same position mutation, reset MutLoc to 0 to "forget"
        %last mutation location.
        if AllowSameMutPos == 1 
            MutLoc = MutLoc*0;
        end
        
        MutCount = 1;
        while MutCount <= SHMrate
            %Find the A, C, G, T locations that are left over for mutations
            SeqRemLoc{1} = (Seq == 'A') & (MutLoc == 0);
            SeqRemLoc{2} = (Seq == 'C') & (MutLoc == 0);
            SeqRemLoc{3} = (Seq == 'G') & (MutLoc == 0);
            SeqRemLoc{4} = (Seq == 'T') & (MutLoc == 0);
            
            %Calculation probability to mutate a nt composition
            BaseCt = [sum(SeqRemLoc{1}) sum(SeqRemLoc{2}) sum(SeqRemLoc{3}) sum(SeqRemLoc{4})];
            BasePropensity = cumsum(BaseCt .* Pnt);
            BasePropensity = BasePropensity / max(BasePropensity);

            %Select the nt composition to mutate.
            Rcomp = find(rand(1) <= BasePropensity);
            Rcomp = Rcomp(1);

            %Determine the new nt
            RnewComp = find(rand(1) <= CPnt2nt(:,Rcomp));
            RnewComp = RnewComp(1);
            
            %Select the nt position of the composition to mutate.
            SeqRemIdx = find(SeqRemLoc{Rcomp} == 1);
            Rpos = SeqRemIdx(ceil(length(SeqRemIdx)*rand(1)));

            %Make sure it does not create a stop codon
            SeqT = Seq;
            SeqT(Rpos) = int2nt(RnewComp);
            AA = nt2aa(SeqT,'ACGTonly','false','frame',ReadFrame);
            if ~isempty(regexp(AA,'*','once'))
                continue
            elseif Rpos >= CDR3start && Rpos <= CDR3start+2; %make sure it doesn't get rid of 104C
                Caa = nt2aa(SeqT(CDR3start:CDR3start+2));
                if ~strcmpi(Caa,'C')
                    continue
                end
            elseif Rpos >= CDR3end-2 && Rpos <= CDR3end; %make sure it doesn't get rid of 118W/F
                Waa = nt2aa(SeqT(CDR3end-2:CDR3end));
                if ~(strcmpi(Waa,'W') || strcmpi(Waa,'F'))
                    continue
                end                
            else
                Seq = SeqT;
                MutLoc(Rpos) = 1;
                MutCount = MutCount + 1;
            end
        end
        
        %Update the new entry after accumulating enough mutations
        NewVDJdata(q,:) = VDJdata(j,:);
        NewVDJdata{q,SeqLoc} = Seq;
        
        %Saving the ancestral seq
        if g == 1
            NewVDJdata{q,RefSeqLoc} = VDJdata{j,SeqLoc}; %Germline sequence
        else
            NewVDJdata{q,RefSeqLoc} = NewVDJdata{q-1,SeqLoc}; %Parent seq
        end
        NewVDJdata{q,SeqNumLoc} = SeqNum + size(VDJdata,1)*g;
        NewVDJdata{q,SeqNameLoc} = SeqNum + size(VDJdata,1)*g;
        q = q + 1;
    end
end
VDJdata  = NewVDJdata;   

%Fill in the details now
VDJdata = buildVDJalignment(VDJdata,NewHeader); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 


DotLoc = find(FileName == '.');
SaveName = sprintf('%s_GrpSize%0.0f_SHMrate%0.0f',FileName(1:DotLoc(end)-1),GrpSize,SHMrate);

%Restructure the alignment data
VDJdata = reformatAlignment(VDJdata,1);
if ispc
    xlswrite([FilePath SaveName '.xlsx'],cat(1,NewHeader,VDJdata));
    xlswrite([FilePath SaveName '.xlsx'],Pnt2nt,'SHMmatrix');
else
    writeDlmFile(cat(1,NewHeader,VDJdata),[FilePath SaveName '.csv'],'\t');
    writeDlmFile(Pnt2nt,[FilePath SaveName 'SHMmatrix.csv'],'\t');
end