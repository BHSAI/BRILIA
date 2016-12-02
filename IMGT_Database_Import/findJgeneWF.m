%findJgeneFW will look for the 118 ForW residue within the J gene reference
%sequences, and return the number of NTs counting from the 5' end to the
%1st letter of the codon that codes for the C residue.
%
%Note: This only has to be done once when setting up the
%IMGT_MOUSE_VDJgene_AllStrain.mat file using the "generateIMGTmat.m"
%script.
%
%  [Cloc, RF, AA] = findJgeneWF(Seq)
%    Cloc = position from the 3' end. 
%      EX: Seq(Cloc:Cloc+2) = codon for W or F.
%
%The C distance from the end must be <= 8. 104C of V, 105, 106, 107,
%108, 109, 110, 111 .... 112, 113, 114, 115, 116, 117, 118W of J

function [Cloc, RF, AA] = findJgeneWF(Seq)
ConservedSeq1 = 'FGXGTT';    
ConservedSeq2 = 'WGXGTT';    

Cloc = zeros(size(Seq,1),1);
RF = zeros(size(Seq,1),1);
AA = cell(size(Seq,1),1);

for j = 1:size(Seq)    
    TempNT = Seq{j};
    TempAA = nt2aa(TempNT,'ACGTonly',0,'frame','all');
    ScoreMat = zeros(3,4);
    for k = 1:3
        [Score1, ~, StartAt1] = convolveSeq(ConservedSeq1,TempAA{k},0,0);
        [Score2, ~, StartAt2] = convolveSeq(ConservedSeq2,TempAA{k},0,0);
        
        %Determine if it's a F or W;
        if Score1(3) > Score2(3) %It's an F
            BestScore = Score1(3);
            BestStartAt = StartAt1(2);
        else %It's an W or tied
            BestScore = Score2(3);
            BestStartAt = StartAt2(2);            
        end
        
        %Make sure * does not occur AFTER the StartAt position
        StopLoc = find(TempAA{k} == '*');
        if isempty(StopLoc); StopLoc = 0; end
        
        ScoreMat(k,1) = BestScore;
        ScoreMat(k,2) = BestStartAt;
        ScoreMat(k,3) = StopLoc(end);
        ScoreMat(k,4) = k;
    end
    
    %Delete those with stop codon after the StartAt position
    DelThis = ScoreMat(:,3) > ScoreMat(:,2);
    ScoreMat(DelThis,:) = [];
    
    BestLoc = find(ScoreMat(:,1) == max(ScoreMat(:,1)));
    BestRF = ScoreMat(BestLoc(1),end);
  
    %Calculate the NT position where the FW is coded
    Cloc(j) = (BestRF-1) + 3*(ScoreMat(BestLoc,2)-1) + 1;
    RF(j) = BestRF;
    AA{j} = TempAA{BestRF};
end

