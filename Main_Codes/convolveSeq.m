%convolveSeq will take only NT sequence 1 and 2, match them,
%and return Score, Alignment, StartAt values. 
%
%  [Score, Alignment, StartAt, MatchAt] = convolveSeq(Seq1, Seq2,
%  [OverhangDir], [TrimEnds])  where Seq1 <= Seq2 in length. Seq2 is
%  normally the longer reference sequence.
%    - Score = [(# of matches); (length of matching nts); (FitScore)];
%        FitScore = sum((cont matches)^2) - sum((cont mismatches)^2).
%    - Alignment = 3xN character array of Seq1 and Seq2 match. 
%    - StartAt = [Seq1StartLoc; Seq2StartLoc], where Seq1StartLoc = 1 by
%        default. Seq2StartLoc can be negative, indicating how many NTs
%        overhang Seq1 over Seq2.
%    - OverhangDir is optional. Forces a preferred Seq1 overhang.
%    - TrimThis is 1 by default. Trims poorly aligned ends.
%
%  convolveSeq(Seq1,Seq2,OverhangDir) will assume Seq1 has a
%  left overhang (OverhangDir = -1), right overhang (OverhangDir = 1), or
%  any overhang (OverhangDir = 0); *If OverhangDir = 'match', then it does
%  NOT do any overhang match, and does a direct matching of Seq1 to Seq2,
%  assuming they have the same lengths.
%
%  Example of Trim for poorly matched ends:
%    Seq1 = 'TAATAATTAAT'
%    Seq2 = 'TCCTAATTGGT'
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0,1)
%    Alignment =
%         TAATAATTAAT
%            |||||   
%         TCCTAATTGGT
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0,0)
%    Alignment =
%         TAATAATTAAT
%         |  |||||  |
%         TCCTAATTGGT
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,1,1)
%         TAATAATTAAT
%         |  |||||   
%         TCCTAATTGGT
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,-1,1)
%         TAATAATTAAT
%            |||||  |
%         TCCTAATTGGT

%  Example of OverhandDir effect:
%    Seq1 = 'CTTAGGAA'
%    Seq2 = 'GGAACTTAGGAACTTA'
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,-1)
%    Alignment =
%         CTTAGGAA------------
%             ||||            
%         ----GGAACTTAGGAACTTA
% 
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0)
%         Alignment =
%         ----CTTAGGAA----
%             ||||||||    
%         GGAACTTAGGAACTTA
%
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,+1)
%         Alignment =
%         ------------CTTAGGAA
%                     ||||    
%         GGAACTTAGGAACTTA----
%
%  Example of Match effect:
%    Seq1 = 'ACGGT'
%    Seq2 = 'CGGTT'
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0)
%        Alignment =
%        ACGGT-
%         |||| 
%        -CGGTT
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,'match')
%         ACGGT
%           | |
%         CGGTT
%
%  Tie breakers:
%  1) If there are equal match/mismatch numbers, break ties based on
%     longest consecutive match.
%    Seq1 = 'CGGGC'
%    Seq2 = 'CGTGCAAACGGTC'
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0)
%         Alignment =
%         --------CGGGC
%                 ||| |
%         CGTGCAAACGGTC

%  2) If are still ties, win based on overhang direction.
%    Seq1 = 'CCCC'
%    Seq2 = 'CCCCTTTTCCCCTTTTCCCC'
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,-1)
%         Alignment =
%         CCCC----------------
%         ||||                
%         CCCCTTTTCCCCTTTTCCCC
%
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,+1)
%         Alignment =
%         ----------------CCCC
%                         ||||
%         CCCCTTTTCCCCTTTTCCCC
%
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0)
%         Alignment =
%         --------CCCC--------
%                 ||||        
%         CCCCTTTTCCCCTTTTCCCC
%
%    X matching allowed here for nucleutides.
%    Seq1 = 'AGXTXXT';
%    Seq2 = 'CCCAGTAGTCCC';
%    [Score, Alignment, StartAt] = convolveSeq(Seq1,Seq2,0)
%         ---AGTXXT---
%            ||||||   
%         CCCAGTAGTCCC
%  * If there are still ties, only select first occurence.
%
%    Seq1 = '******GACTGG';
%    Seq2 = 'GACTATGACTGG';
%    
%
%  Allow special scoring for mutations
%    Seq1 = 'ACCGAT';
%    Seq2 = 'ACTGGT';
%    [Score, ~, ~] = convolveSeq(Seq1,Seq2,0,1,0);
%        Score(3) = 4;
%    [Score, ~, ~] = convolveSeq(Seq1,Seq2,0,1,1);
%        Score(3) = 9;
%
function varargout = convolveSeq(Seq1,Seq2,varargin)
%Establish default alignment settings
OverhangDir = 0; %Default 0, best match without any overhang preference. 
TrimEnds = 1; %Default 1, yes trim, but the side changes with OverhangDir. EX: If OverhangDir = 0 , Trim BOTH ends. OverhangDir = 1, trim right only. OverhangDir = -1, trim left only.
AllowedMiss = 0; %This affects the calcAlignScore

%All the options possible
OverhangModes = {'LeftPad' 'Left' 'Any' 'Right' 'RightPad' 'Match' 'LeftPad'} ;
OverhangModeNum = [-3 -1 0 1 -3 2];
TrimEndsModes = {'none' 'trim' };
TrimEndsModeNum = [0 1];

%Modify defaults based on inputs. NOTE: Inputting number mode is faster
%than text modes!
if length(varargin) >= 1
    if ischar(varargin{1})
        for k = 1:length(OverhangModes)
            if strcmpi(OverhangModes{k},varargin{1})
                OverhangDir = OverhangModeNum(k);
                break
            end
        end
    else
        OverhangDir = varargin{1};
    end
    
    if length(varargin) >= 2        
        if ischar(varargin{2})
            for k = 1:length(TrimEndsModes)
                if strcmpi(TrimEndsModes{k},varargin{2})
                    TrimEnds = TrimEndsModeNum(k);
                    break
                end
            end
        else
            TrimEnds = varargin{2};
        end
    
        if length(varargin) >= 3
            AllowedMiss = varargin{3}; %For the Leneient score mode.
        end
    end
end

%Error checking
if OverhangDir == 2
    if length(Seq1) ~= length(Seq2)
        error('For option ''match'', length of Seq1 and Seq2 must be same.');
    end
end
if iscell(Seq1) 
    error('Seq1 should be a char array');
end

%Format input to standardized format
Seq1 = upper(Seq1);
Seq2 = upper(Seq2);
SeqMap = 1:(size(Seq1,2)+size(Seq2,2)); %Used to convert binary search to numeric indexes.

%Check for existence of X in Seq2 anywhere;
HaveX = 0;
Xcheck2 = regexp(Seq2,'X','once');
if ~isempty(Xcheck2)
    HaveX = 1;
end
if HaveX == 0 %Do a final check on Seq1
    for q = 1:size(Seq1,1)
        Xcheck1 = regexp(Seq1(1,:),'X','once');
        if ~isempty(Xcheck1)
            HaveX = 1;
            break
        end
    end
end

%Determine the shorter sequences & adjust OverhandDir too. Code works best
%if Seq1 is shorter t han Seq2
if size(Seq1,2) <= size(Seq2,2)
    SeqB = Seq2;
    FlipThis = 0;
else
    SeqA = Seq2;
    OverhangDir = -OverhangDir;
    FlipThis = 1; %Going to have to flip the results at the end.
    %Also scale down the AllowedMiss for the short Seq
    AllowedMiss = ceil(AllowedMiss / length(Seq1) * length(Seq2));
end

%==========================================================================
%Enabling multiple alignment of trimmed, Seq1
AllScores = cell(size(Seq1,1),1);
AllAlignments = cell(size(Seq1,1),1);
AllStartAts = cell(size(Seq1,1),1);
AllMatchAts = cell(size(Seq1,1),1);

for q = 1:size(Seq1,1)
    if FlipThis == 0
        SeqA = Seq1(q,:);
    else
        SeqB = Seq1(q,:);
    end

    %Determine if there are any dummy * characters in any sequenc.
    StarA = sum((SeqA == '*'));
    StarB = sum((SeqB == '*'));
    MaxOvlp = min([(length(SeqA) - StarA) (length(SeqB) - StarB)]);
    
    if q == 1 %You only do this once
        if OverhangDir == 2 
            SeqT = SeqB;
            s1 = length(SeqA);
        else

            %Determine the search matrix (avoids the for loop checking)
            if OverhangDir == -1
                s1 = 1;
                s2 = length(SeqA); %It's SeqA because otherwise, you don't force overhang matching.
            elseif OverhangDir == 1
                s1 = length(SeqB);
                s2 = length(SeqA) + length(SeqB) - 1;
            else
                s1 = 1;
                s2 = length(SeqA) + length(SeqB) - 1;
            end
            SeqT = makeDiagonalSeq(SeqA,SeqB,s1,s2);
        end

        if HaveX == 1
            SeqTn = (SeqT == 'X'); %This is used later to overrule "X" nt mismatch. X should match with any A G C T.
        end
        
        %Get the number of NTs being matched per row of SeqT
        SeqOvlp = sum((SeqT ~= '-'),2);
        SeqOvlp(SeqOvlp > MaxOvlp) = MaxOvlp;
    end
    
    SeqM = repmat(SeqA,size(SeqT,1),1);
    
    if HaveX == 1
        SeqMn = (SeqM == 'X'); %This is used to overrule any "X" nt mismatch. X should match with ant A G C T.
        SeqTM = (SeqT == SeqM) | SeqMn | SeqTn; %Matched Results, allowing for X matching too.
        SeqDash = SeqT == '-'; %Gets rid of matching the X to '-'.
        SeqTM(SeqDash) = 0;
    else
        SeqTM = (SeqT == SeqM);
    end
    
    SumSeqTM = sum(SeqTM,2);
    SeqIdentity = SumSeqTM.^2-(SeqOvlp-SumSeqTM).^2;  %Should give maximum possible score.
    [SeqIdentity, SeqIdentityIdx] = sort(SeqIdentity,'descend');
    
    %Find the top 5 unique identities, avoiding using unique function
    UnqIdentityScore = SeqIdentity(1);
    UnqIdentityNum = 1;
    TopMatchNum = 5;
    for w = 2:length(SeqIdentity);
        if UnqIdentityScore ~= SeqIdentity(w)
            UnqIdentityNum = UnqIdentityNum + 1;
            UnqIdentityScore = SeqIdentity(w);
        end
        if UnqIdentityNum > TopMatchNum
            break
        end
        if SeqIdentity(w) < 0; break; end
    end
    EvalThese = SeqIdentityIdx(1:w);
    if isempty(EvalThese); EvalThese = 1; end %caused of OverhangDir == 2.
    
    %Determine which 1st and last match occurrence
    FitMatrix = zeros(size(SeqTM,1),3);
    for b = 1:size(EvalThese,1)
        j = EvalThese(b);
        MatchIdx = SeqMap(SeqTM(j,:));
        if isempty(MatchIdx); continue; end
        
        %Trim bad ends once EX: |__|__||||____| will become ||||
        if TrimEnds == 1
            if length(MatchIdx) >= 3;
                if OverhangDir == -1 %Trim left ends
                    MatchedNts = double(SeqTM(j,MatchIdx(1):end));
                    MatchedNts(MatchedNts == 0) = -1;
                    SumMatch = cumsum(MatchedNts);
                    SumMatch(MatchedNts > 0) = 1; %Ensure the cut loc is at a mismatches
                    CutLoc = find(SumMatch <= 0);
                    if ~isempty(CutLoc)
                        CutLoc = CutLoc(end);
                        DelThese = MatchIdx(MatchIdx <= MatchIdx(1)+CutLoc-1);
                        SeqTM(j,DelThese) = 0;
                        MatchIdx(MatchIdx <= MatchIdx(1)+CutLoc-1) = [];
                    end
                                
                elseif OverhangDir == 1 %Trim right ends
                    MatchedNts = double(SeqTM(j,1:MatchIdx(end)));
                    MatchedNts(MatchedNts == 0) = -1;
                    MatchedNts = fliplr(MatchedNts);
                    SumMatch = cumsum(MatchedNts);
                    SumMatch(MatchedNts > 0) = 1;
                    CutLoc = find(SumMatch <= 0);
                    if ~isempty(CutLoc)
                        CutLoc = CutLoc(end);
                        DelThese = MatchIdx(MatchIdx >= MatchIdx(end)-CutLoc+1);
                        SeqTM(j,DelThese) = 0;
                        MatchIdx(MatchIdx >= MatchIdx(end)-CutLoc+1) = [];
                    end
                               
                else
                    if (MatchIdx(2) - MatchIdx(1)) >= 3 %Contains a double gap
                        SeqTM(j,MatchIdx(1)) = 0; %1st, set to 0. This is used later for alignment.
                        MatchIdx(1) = MatchIdx(2); %2nd, temp change for FitScore calc
                    end
                    if (MatchIdx(end) - MatchIdx(end-1)) >= 3 %Contains a double gap
                        SeqTM(j,MatchIdx(end)) = 0; %1st, set to 0. This is used later for alignment.
                        MatchIdx(end) = MatchIdx(end-1); %2nd, temp change for FitScore calc
                    end
                end
            end
        end
        
        if isempty(MatchIdx); continue; end
        
        FitMatrix(j,2) = MatchIdx(end) - MatchIdx(1) + 1;
        
        %Make the alignment 3xN char array early
        if OverhangDir == -1
            g1 = MatchIdx(1);
            g2 = size(SeqTM,2);
        elseif OverhangDir == 1
            g1 = 1;
            g2 = MatchIdx(end);
        else
            g1 = MatchIdx(1);
            g2 = MatchIdx(end);            
        end
        MatchedResults = SeqTM(j,g1:g2);
            
%         TopSeq = SeqA(MatchIdx(1):MatchIdx(end));
%         BotSeq = SeqT(j,MatchIdx(1):MatchIdx(end));
%         MidSeq = repmat(' ',1,length(TopSeq));
%         MidSeq(SeqTM(j,MatchIdx(1):MatchIdx(end))) = '|';
%         if FlipThis == 0
%             Alignment = [TopSeq; MidSeq; BotSeq];
%         else
%             Alignment = [BotSeq; MidSeq; TopSeq];
%         end
%         
        FitMatrix(j,3) = calcAlignScore(MatchedResults,AllowedMiss);
    end
    
    FitMatrix(:,1) = sum(SeqTM,2); %Calculate sum of matches
    BestMatch = find(FitMatrix(:,3) == max(FitMatrix(:,3)));

    %If there's a tie, win base on OverhangDir
    if length(BestMatch) >  1
        if OverhangDir == 1 %Take one with more left overhang
            BestMatch = BestMatch(end);
        elseif OverhangDir == -1 %Take one with more right overhang
            BestMatch = BestMatch(1);
        else %Look for "middle" BestMatch
            BestMatch = BestMatch(round(length(BestMatch)/2));
        end
    end

    %Determine the start loc of Seq2, and final alignment length
    Seq2start = BestMatch - (length(SeqA) - 1) + (s1 - 1); %Negative values possible.
    if Seq2start <= 0; Seq2start = Seq2start - 1; end %This is done since there is no "0" start loc.
    if Seq2start < 0;
        TotalLen = abs(Seq2start) + length(SeqB);
    else
        TotalLen = Seq2start + length(SeqA)-1;
        if TotalLen < length(SeqB)
            TotalLen = length(SeqB);
        end       
    end

    AllScores{q} = FitMatrix(BestMatch,1:3)';
    
    if nargout >= 2
        %Make the alignment 3xN char array
        Alignment = repmat(['-';' ';'-'],1,TotalLen); 
        MatchPat(1:length(SeqA)) = ' ';
        MatchPat(SeqTM(BestMatch,:)) = '|';
        if Seq2start < 0
            Alignment(1,1:length(SeqA)) = SeqA;
            Alignment(2,1:length(SeqA)) = MatchPat;
            Alignment(3,abs(Seq2start)+1:abs(Seq2start)+length(SeqB)) = SeqB;
        else
            Alignment(1,Seq2start:Seq2start+length(SeqA)-1) = SeqA;
            Alignment(2,Seq2start:Seq2start+length(SeqA)-1) = MatchPat;
            Alignment(3,1:length(SeqB)) = SeqB;
        end

        %Flip the stuff
        if FlipThis == 1
            Alignment = flipud(Alignment);
            Seq2start = -Seq2start + 1;
            if Seq2start == 0
                Seq2start = 1;
            end
        end

        %format output
        AllAlignments{q} = Alignment;
        AllStartAts{q} = [1; Seq2start];
        if nargout == 4
            MatchIdx = SeqMap(Alignment(2,:)=='|');
            if isempty(MatchIdx)
                MatchIdx = [0 0];
            end
            AllMatchAts{q} = [MatchIdx(1); MatchIdx(end)];
        end
    end
end

%Format outputs
if size(Seq1,1) == 1
    if nargout >= 1
        varargout{1} = AllScores{1};
        if nargout >= 2
            varargout{2} = AllAlignments{1};
            if nargout >= 3
                varargout{3} = AllStartAts{1};
                if nargout >= 4
                    varargout{4} = AllMatchAts{1};
                end
            end
        end
    end
else %Save outputs as cells
    if nargout >= 1
    varargout{1} = AllScores;
        if nargout >= 2
            varargout{2} = AllAlignments;
            if nargout >= 3
                varargout{3} = AllStartAts;
                if nargout >= 4
                    varargout{4} = AllMatchAts;
                end
            end
        end
    end
end
