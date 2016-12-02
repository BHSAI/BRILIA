%makeClassifier takes a Seq, VMDNJ lengths, V and D and J alignment results
%to make the classifier. Classifier identifies which NT in the seq belongs
%to which gene, and whether or not it matches with the reference gene.
%
%  VDJdata = makeClassifier(VDJdata,NewHeader)
%    where Seq = NT sequence; VMDNJ is a 1x5 matrix of V, N2, D, N1, J
%    lengths; V/D/J align are the alignment results between sample and V D
%    or J genes. 
%
%    Example Classifier = 'VVVVVVVVVVVVVVVVppmmmppDDDDDDppnnppJJJJJJJJJ';
%    V = matched with V reference
%    v = within V matching region, mismatched though
%    p = p-nucleotides between the VD region
%    m = VD junction N2 nucleotides
%    D = matched with D reference
%    d = within D matching region, mismatched though
%    b = p-nucleotides between the DJ region
%    n = DJ junction N1 nucleotides
%    J = matched with J reference
%    j = within J matching region, mismatched though

function VDJdata = makeClassifier(VDJdata,NewHeader)
getHeaderVar;

%Correct the classifiers according to the group, mainly p and b's only.
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
GrpNumUnq = unique(GrpNum);
for y = 1:length(GrpNumUnq)
    IdxLoc =  find(GrpNumUnq(y) == GrpNum);
    j = IdxLoc(1);
    
    %Extract information
    Seq = VDJdata{j,SeqLoc};
    RefSeq = VDJdata{j,RefSeqLoc};
    VMDNJ = cell2mat(VDJdata(j,LengthLoc));
    
    if sum(VMDNJ) ~= length(Seq) || sum(VMDNJ) ~= length(RefSeq)
        continue;
    end
        
    %Find V classifier
    if VMDNJ(1) > 0
        ClassV = repmat('V',1,VMDNJ(1));
    else
        ClassV = '';
    end
    
    %Find N2 classifier + p nucleotides
    if VMDNJ(2) > 0
        ClassM = repmat('M',1,VMDNJ(2));
    else
        ClassM = '';
    end
    
    %Find D classifier
    if VMDNJ(3) > 0
        ClassD = repmat('D',1,VMDNJ(3));
    else
        ClassD = '';
    end
        
    %Find N1 classifier
    if VMDNJ(4) > 0
        ClassN = repmat('N',1,VMDNJ(4));
    else
        ClassN = '';
    end

    %Find J classifier
    if VMDNJ(5) > 0
        ClassJ = repmat('J',1,VMDNJ(5));
    else
        ClassJ = '';
    end
    
    %Find the VD and DJ p nucleotides, marked as "p" and "b" respectively.
    Vref3del = VDJdata{j,DelLoc(1)};
    Dref5del = VDJdata{j,DelLoc(2)};
    Dref3del = VDJdata{j,DelLoc(3)};
    Jref5del = VDJdata{j,DelLoc(4)};
    
    if Vref3del == 0 && VMDNJ(2) > 0
        for k = 1:VMDNJ(2)
            s1 = VMDNJ(1) - k + 1;
            s2 = VMDNJ(1) + k;
            if s1 < 1; continue; end
            if Seq(s1) ~= seqcomplement(Seq(s2));
                break
            else
                ClassM(k) = 'B';
            end
        end
    end
    
    if Dref5del == 0 && VMDNJ(2) > 0
        for k = 1:VMDNJ(2)
            s1 = sum(VMDNJ(1:2)) - k + 1;
            s2 = sum(VMDNJ(1:2)) + k;
            if Seq(s1) ~= seqcomplement(Seq(s2));
                break
            else
                ClassM(end-k+1) = 'B';
            end
        end
    end
     
    if Dref3del == 0 && VMDNJ(4) > 0
        for k = 1:VMDNJ(4)
            s1 = sum(VMDNJ(1:3)) - k + 1;
            s2 = sum(VMDNJ(1:3)) + k;
            if Seq(s1) ~= seqcomplement(Seq(s2));
                break
            else
                ClassN(k) = 'P';
            end
        end
    end   
    
    if Jref5del == 0 && VMDNJ(4) > 0
        for k = 1:VMDNJ(4)
            s1 = sum(VMDNJ(1:4)) - k + 1;
            s2 = sum(VMDNJ(1:4)) + k;
            if s2 > sum(VMDNJ); continue; end
            if Seq(s1) ~= seqcomplement(Seq(s2));
                break
            else
                ClassN(end-k+1) = 'P';
            end
        end
    end   
    
    AllClass = sprintf('%s%s%s%s%s',ClassV,ClassM,ClassD,ClassN,ClassJ);
    
    %Fillin the Classifier and Formatted Seq 
    for k = 1:length(IdxLoc)
        CurSeq = VDJdata{IdxLoc(k),SeqLoc};
        MissLoc = CurSeq ~= RefSeq;
        CurClass = AllClass;
        CurClass(MissLoc) = lower(CurClass(MissLoc));
        
        VDJdata{IdxLoc(k),FormClassLoc(1)} = formatSeq(VDJdata{IdxLoc(k),SeqLoc},CurClass);
        VDJdata{IdxLoc(k),FormClassLoc(2)} = CurClass;
    end
    clear ClassV ClassM ClassD ClassN ClassJ
end

