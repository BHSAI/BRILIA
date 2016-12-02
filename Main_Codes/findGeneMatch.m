%Perform gene matching using Xmap database on NTsample, looking for similar
%nucleotide residues. This is the core alignment-based, gene search.
%
%  AllResults = findGeneMatchNT(Seq,Xmap,VDJ) where NTsample is a
%  trimmed, character or cell array of sequence, Xmap is a cell matrix
%  where column 1 contains the refence sequences, VDJ is either "V" or "D"
%  or "J" to denote the gene category. 
%  
%  AllResults = findGeneMatchNT(Seq,Xmap,VDJ,['normal' or 'shm']) will
%  change the scoring mode, such that 'shm' mutation mismatches does not
%  result in too low of a alignment score. Useful for annotating D's. (See
%  calcAlignScore for details on how this affects scoring)
%
%  AllResults is a cell array of output of size(Seq1,1) x 6. Each Col
%  contains:
%    Col1: Gene index, or row, number in Xmap
%    Col2: Full gene name
%    Col3 = [L M R] of the reference gene
%    Col4 = [L M R] of the Seq
%    Col5 = Identity match [MatchNum TotMatchPossible]
%    Col6 = Alignment text 

function AllResults = findGeneMatch(Seq,Xmap,VDJ,varargin)
%Make sure NTsample is not a cell
if iscell(Seq)
    Seq = char(Seq);
end

%Setup the overhang directions depending on the gene type.
if strcmpi(VDJ,'V')
    OverhangDir = 1;
    TrimEnds = 1;
elseif strcmpi(VDJ,'J')
    OverhangDir = -1;
    TrimEnds = 1;
elseif strcmpi(VDJ,'D')
    OverhangDir = 0;
    TrimEnds = 1; 
end

%Determine AllowedMiss in alignment
if ~isempty(varargin)
    AllowedMiss = varargin{1};
else
    AllowedMiss = 0;
end

%Create a score for each alignment. %Score NThit NTtot SampLMR GeneLMR
%LMR = [LeftflankNTs MidNTs RightflankNTs] by lengths. EX: [1 3 1];

%Determine best fit via NT matching
FitScore = zeros(size(Seq,1),size(Xmap,1)); %cat(2,zeros(size(Xmap,1),1),[1:size(Xmap,1)]');
for w = 1:size(Xmap,1)
    %Determine best fit via NT matching
    if isempty(Xmap{w,1}); continue; end
    FitScoreT = convolveSeq(Seq,Xmap{w,1},OverhangDir,TrimEnds,AllowedMiss);
    if iscell(FitScoreT);
        FitScoreT = cell2mat(FitScoreT); %This will convert a cell in a cell into a linear vector,  hence must extract every 3rd to get score.
    end
    FitScore(:,w) = FitScoreT(3:3:end);
end

%Determine the column with max score, for each seq 
MaxScores = max(FitScore,[],2);
GeneIdx = cell(size(Seq,1),1); %This is an output
for j = 1:size(Seq,1)
    GeneIdx{j} = find(MaxScores(j) == FitScore(j,:)); 
end

%Redo the alignments to extract the LMR + alignment infos
AllResults = cell(size(Seq,1),6);
for j = 1:size(Seq,1)
    BestMatch = GeneIdx{j};
    NumMatch = length(BestMatch);
    AlignResults = cell(NumMatch,6);
    for k = 1:length(BestMatch)
        w = BestMatch(k);
        [AlignScore, Alignment, ~, MatchAt] = convolveSeq(Seq(j,:),Xmap{w,1},OverhangDir,TrimEnds,AllowedMiss); 
        NThit = sum(Alignment(2,:) == '|');

        %Extract the relevant informations of alignment
        SampSeq = Alignment(1,:);
        SampLoc = find(SampSeq ~= '-');
        SampLoc = SampLoc(1,[1 end]);

        GeneSeq = Alignment(3,:);
        GeneLoc = find(GeneSeq ~= '-');
        GeneLoc = GeneLoc(1,[1 end]);

        %Determine the LMRs flanking lengths
        Ls = MatchAt(1) - SampLoc(1);
        Ms = MatchAt(end) - MatchAt(1) + 1;
        Rs = SampLoc(end) - MatchAt(end);
        Lg = MatchAt(1) - GeneLoc(1);
        Mg = Ms;
        Rg = GeneLoc(end) - MatchAt(end);   

        %Fill in AlignResults
        AlignResults{k,1} = w; %GeneNum in the Xmap
        AlignResults{k,2} = Xmap{w,3}; %GeneName
        AlignResults{k,3} = [Lg Mg Rg];
        AlignResults{k,4} = [Ls Ms Rs];
        AlignResults{k,5} = [NThit AlignScore(3)];
        AlignResults{k,6} = Alignment;
    end    
    AlignResults = reduceFamily(AlignResults);

    %Now find the match with the MINIMUM total deletion in the appropriate
    %direction
    LMRgene = cell2mat(AlignResults(:,3));
    if VDJ == 'V'
        MinLoc = find(LMRgene(:,3)==min(LMRgene(:,3)));
    elseif VDJ == 'J'
        MinLoc = find(LMRgene(:,1)==min(LMRgene(:,1)));
    else
        Max5and3 = max(LMRgene(:,[1 3]),[],2);
        MinLoc = find(Max5and3 == min(Max5and3));
    end

    %Report only 1st results if there are still ties.
    AllResults(j,:) = AlignResults(MinLoc(1),:);
end