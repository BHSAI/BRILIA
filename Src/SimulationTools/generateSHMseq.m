%generateSHMseq will introduce SHM per each sequence clonal group. If the
%input contains grouped sequences, will randomly select one sequence and
%generate AddSize number of linear descendants. Hence, running
%generateSHMseq on the same VDJdata will create branched lineages.
%
%  VDJdata = generateSHMseq(VDJdata,VDJheader,SHMperc,AddSiz,AllowRepeatMut)
%
%  INPUT
%    SHMperc: % of nts that must undergo SHM in order to save that seq
%    BranchLength: # of seq to add linearly per branch and per group.
%    BranchCount: # of lineage branches to add per group.
%    AllowRepeatMut ['y' or 'n']: Allows SHM to occur in the same location
%       as a previous SHM.
%
%  OUTPUT
%    Tdata: cell matrix containing VDJ seq AND descendant seq added.
%
%  NOTE
%    VDJdata will be expanded by (AddSize * CloneCount) entries, where
%    CloneCount is the number of unique clusters in VDJdata. CloneCount
%    does not have to be specified, as it is determined via unique GrpNum
%    in VDJdata.

function VDJdata = generateSHMseq(VDJdata,VDJheader,SHMperc,BranchLength,BranchCount,AllowRepeatMut)
%Determine chain
H = getHeavyHeaderVar(VDJheader);
L = getLightHeaderVar(VDJheader);
if H.SeqLoc > 0 && L.SeqLoc > 0
    Chain = 'HL';
elseif H.SeqLoc > 0
    Chain = 'H';
else
    Chain = 'L';
end

%Determine if mutation can occur at the same position
AllowRepeatMut = lower(AllowRepeatMut);

%Calc probablity nt goes from (column) nt to (row) nt. 
%From A Collins 2015 data, analyzed by BRILIA Dev 1.9.0.
% %         A0      C0      G0      T0
% Pnt2nt = [0       10395   25057   9981      %A1
%           12809   0       9393    24737     %C1
%           34286   10005   0       7375      %G1
%           15368   30254   5561    0    ];   %T1

%Calc probablity nt goes from (column) nt to (row) nt. 
%From D Lee 2017 BRILIA mouse data, processed by BRILIA v1.13.0
%         A0      C0     G0      T0
Pnt2nt = [0       7320   22093   5798    %A1
          9228    0      8200    19160   %C1
          30650   7874   0       7492    %G1
          12086   23588  5886    0    ]; %T1
Pnt2nt = Pnt2nt./sum(Pnt2nt(:)); 
Pnt = sum(Pnt2nt,1); %Probablity that a nt type is mutated.
CPnt2nt = cumsum(Pnt2nt,1)./repmat(sum(Pnt2nt,1),4,1); %Cumulative probability per each pairwise mutation

%Determine the group numbers
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);

%Create the new size VDJdata
CloneCount = length(UnqGrpNum);
SeqAddCount = CloneCount * (1 + BranchCount * BranchLength);
NewVDJdata = cell(SeqAddCount, size(VDJdata,2));

q = 1; %Counter for NewVDJdata
for y = 1:length(UnqGrpNum) %Find group and select random one to be the ancestor sequence
    %Transfer over current data
    GrpIdx = find(GrpNum == UnqGrpNum(y)); %Find Group
    GrpSize = length(GrpIdx); %For keeping track of group growth
    NewGrpIdx = q:q+GrpSize-1;
    NewVDJdata(NewGrpIdx,:) = VDJdata(GrpIdx,:);
    q = q + length(GrpIdx);

    %For each branch, compute linear lineages
    for b = 1:BranchCount
        SeqRelIdx = ceil(GrpSize * rand(1));
        SeqIdx = NewGrpIdx(1) + SeqRelIdx - 1; %Asign random ancestor sequence in the NewVDJdata seq group

        %For each chain, perform SHMs
        q0 = q; %restart at this position
        for k = 1:length(Chain)
            if strcmpi(Chain(k),'H')
                B = H;
                q = q0;
            else
                B = L;
                q = q0;
            end
                   
            %Extract basics
            Seq = upper(NewVDJdata{SeqIdx,B.SeqLoc});
            SeqNum = NewVDJdata{SeqIdx,B.SeqNumLoc};
            CDR3s = NewVDJdata{SeqIdx,B.CDR3Loc(3)};
            CDR3e = NewVDJdata{SeqIdx,B.CDR3Loc(4)};
            ReadFrame = mod(CDR3e,3)+1;
            SHMcount = round(SHMperc/100 * length(Seq));   

            %Add a mutated seq after complete SHMcount number of SHMs
            MutLoc = zeros(1,length(Seq),'logical'); %Keeps track of where mutations were done
            for g = 1:BranchLength
                %"Forget" MutLoc to allow SHMs in the same positions.
                if AllowRepeatMut == 'y'; 
                    MutLoc = zeros(1,length(Seq),'logical');
                end

                %----------------------------------------------------------
                %Perform mutations 1 by 1 until reaching SHMcount number of SHMs
                MutCount = 0;
                while MutCount < SHMcount
                    %Find the A, C, G, T locations that are allowed to mutate
                    if MutCount == 0 %only do this first time. Have to update SeqRemLoc and BaseCt per iteration!
                        SeqRemLoc = cell(1,4);
                        BaseCt = zeros(1,4);
                        for w = 1:4
                            SeqRemLoc{w} = (Seq == int2nt(w)) & ~MutLoc;
                            BaseCt(w) = sum(SeqRemLoc{w});
                        end
                    end

                    %Calculate probability to mutate a nt composition
                    BasePropensity = cumsum(BaseCt .* Pnt);
                    BasePropensity = BasePropensity / max(BasePropensity);

                    %Select the nt composition to mutate, X0;
                    Rcomp = find(rand(1) <= BasePropensity);
                    Rcomp = Rcomp(1);

                    %Determine the new nt, X1;
                    RnewComp = find(rand(1) <= CPnt2nt(:,Rcomp));
                    RnewComp = RnewComp(1);

                    %Determine which X1 nt to mutation, by position
                    SeqRemIdx = find(SeqRemLoc{Rcomp} == 1);
                    Rpos = SeqRemIdx(ceil(length(SeqRemIdx)*rand(1)));

                    %Test out the new mutation
                    SeqT = Seq;
                    SeqT(Rpos) = int2nt(RnewComp);

                    %Determine absolution position of full codon affected
                    UnusedNT = ReadFrame - 1;
                    Ncodons = floor((Rpos - UnusedNT)/3);
                    RemNT = rem((Rpos - UnusedNT),3);
                    if RemNT == 0
                        Ncodons = Ncodons - 1; 
                    end
                    CodonPos = UnusedNT + Ncodons*3 + 1 : UnusedNT + (1+Ncodons)*3;
                    CodonPos(CodonPos <= 0) = [];
                    CodonPos(CodonPos > length(Seq)) = [];

                    %Make sure it does not create a stop codon
                    if length(CodonPos) == 3
                        EvalCodon = SeqT(CodonPos);
                        NewAA = nt2aa(EvalCodon,'ACGTonly',false,'frame',1,'alternative',false);
                        if NewAA == '*' %Has stop codon
                            continue
                        end

                        %Ensure it's not at the CDR3 locations
                        if ~isempty(intersect(CodonPos, CDR3s))
                            if NewAA ~= 'C'
                                continue;
                            end
                        elseif ~isempty(intersect(CodonPos, CDR3e))
                            if NewAA ~= 'W' || NewAA ~= 'F'
                                continue;
                            end
                        end
                    end

                    %Accept SeqT as Seq if mutation is valid, and update variables
                    Seq = SeqT;
                    MutLoc(Rpos) = 1;
                    MutCount = MutCount + 1;
                    SeqRemLoc{Rcomp}(Rpos) = 0; %Remove original letter 
                    SeqRemLoc{RnewComp}(Rpos) = 1; %Add new letter
                    BaseCt(Rcomp) = BaseCt(Rcomp) - 1; %Decrease nt comp
                    BaseCt(RnewComp) = BaseCt(RnewComp) + 1; %Increase nt comp
                end
                %End of mutated sequence generation
                %----------------------------------------------------------

                %Update the new entry after accumulating enough mutations
                if k == 1
                    NewVDJdata(q,:) = NewVDJdata(SeqIdx,:); %Copy existing info
                end
                NewVDJdata{q,B.SeqLoc} = Seq; %Replace Seq with mutated one
                NewVDJdata{q,B.CDR3Loc(1)} = nt2aa(Seq(CDR3s:CDR3e),'ACGTonly',false,'alternative',false); %Replace Seq with mutated one        

                %Save the ancestral seq
                if g == 1
                    NewVDJdata{q,B.RefSeqLoc} = NewVDJdata{SeqIdx,B.SeqLoc}; %First "Ancestor" seq
                else
                    NewVDJdata{q,B.RefSeqLoc} = NewVDJdata{q-1,B.SeqLoc}; %Immediate Parent seq
                end

                %Assign new number to seq (will renumber later anyways)
                NewVDJdata{q,B.SeqNumLoc} = SeqNum + size(VDJdata,1)*g;
                NewVDJdata{q,B.SeqNameLoc} = num2str(SeqNum + size(VDJdata,1)*g);
                q = q + 1;
            end
        end
        GrpSize = GrpSize + BranchLength;
    end
    VDJdata = countSHM(NewVDJdata,VDJheader); %SHM info on the VMDNJ segments, and convert NewVDJdata to VDJdata for output.
end
