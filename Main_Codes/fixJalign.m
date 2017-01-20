%fixJalign will find another J match IF the CDR3end location is >
%SeqLength, but still has a conserved Trp of Phe residue at the end. This
%error occurs often if there is a trimmed sequence ending with the 118W
%instead of going further.

function VDJdata = fixJalign(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end
getHeaderVar;

GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    RootLoc = IdxLoc(1);
    CDR3start = VDJdata{RootLoc,CDR3Loc(3)};
    CDR3end = VDJdata{RootLoc,CDR3Loc(4)};
    Seq = VDJdata{RootLoc,SeqLoc};
    
    try %No error handling. Just try.
        if CDR3end > length(Seq) %Check if there's something wrong with the J alignment
            CDR3seq = Seq(CDR3start:end);
            if mod(length(CDR3seq),3) == 0 %Still a chance it's a bad J alignment
                LastAA = nt2aa(Seq(end-2:end));
                if (LastAA ~= '*') && (strcmpi(LastAA,'W') || strcmpi(LastAA,'F')) %Still a chance it's a bad J alignment
                    Tdata = VDJdata(IdxLoc,:);
                    Jlen = Tdata{1,LengthLoc(5)}; %Which includes the 3 codons for 118residue
                    Jseq = Seq(end-Jlen+1:end);
                    JtempDel = zeros(size(Jmap,1),1);
                    JmapT = Jmap;
                    for k = 1:size(JmapT,1)
                        if isempty(JmapT{k,1}); continue; end
                        WlocEnd = JmapT{k,end}+2; %Locations of the Wloc end codon nt
                        if WlocEnd > Jlen
                           JtempDel(k) = WlocEnd - Jlen;
                           JmapT{k,1} = JmapT{k,1}(JtempDel(k)+1:JtempDel(k)+Jlen);
                           JmapT{k,end} = JmapT{k,end}-JtempDel(k);
                        end
                    end

                    %For all J's, find the one with highest align score.
                    Jmatch = cell(size(JmapT,1),6);
                    DelThis = zeros(size(Jmatch,1),1) > 1;
                    for w = 1:size(Jmatch,1)
                        if isempty(JmapT{w,1}); 
                            DelThis(w) = 1;
                            continue; 
                        end
                        Jmatch(w,:) = findGeneMatch(Jseq,JmapT(w,:),'J');
                        Jmatch{w,1} = w;
                        Jmatch{w,3}(1) = Jmatch{w,3}(1) + JtempDel(w);
                    end
                    Jmatch(DelThis,:) = [];

                    %Find the one with the smallest deletions
                    LMRj = cell2mat(Jmatch(:,3));
                    BestLoc = LMRj(:,1) == min(LMRj(:,1));
                    Jmatch = Jmatch(BestLoc,:);
                    Jmatch = reduceFamily(Jmatch);
                    Jmatch = Jmatch(1,:); %If there's still ties, pick 1.                
                    Tdata(:,[FamNumLoc(3) FamLoc(3) DelLoc(4)]) = repmat([Jmatch(1) Jmatch(2) Jmatch{1,3}(1)],size(Tdata,1),1);

                    %See if it fixes the CDR3 problem
                    Tdata = findCDR3(Tdata,NewHeader); %Get the CDR3 seq and info                
                    if Tdata{1,CDR3Loc(4)} <= length(Seq)
                        disp('correct cdr3')
                        Tdata = buildRefSeq(Tdata,NewHeader,'germline','single'); %must do singles, since group therapy not done.
                        Tdata = buildVDJalignment(Tdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
                        Tdata = makeClassifier(Tdata,NewHeader); %Classifier + FormattedSeq
                        Tdata = appendMutCt(Tdata,NewHeader); %SHM infor on the VMDNJ segments
                        VDJdata(IdxLoc,:) = Tdata;
                    end
                end
            end
        end
    catch
        WarningMsg = sprintf('Warning at %s, sequence group # %d',mfilename,UnqGrpNum(y));
        disp(WarningMsg);
    end
end
