%divideNregion will analyze a sequence in the N region, and determine if it
%was made from TDT via 2-directions or 1-direction.
%
%Note: The assumption is that 3' additions will add G's and A's. If it's on
%the leading strand, you'll see G's and A's, where as if it's on the other
%strand, you'll get C's and T's. 
%
%  Examples:
%  Vague scenarios (no splitting can be decided. Need more information.)
%     Nseq = 'cg';   
%     Nseq = 'gc';    
%     Nseq = 'cgcg';  
%     Nseq = 'gcgc';  
%     Nseq = 'cgcgc'; 
%     Nseq = 'gcgcg'; 
%     Nseq = 'ccccgggg'; It's in the wrong other, as g/a should be on left.
% 
%  Obvious scenarious
%     Nseq = 'gg';     %g | ''
%     Nseq = 'cc';     %''| c
%     Nseq = 'ccgcg';   %''   | ccgc
%     Nseq = 'gcggc';   %gcgg | ''
%     Nseq = 'gggccc'; %ggg | ccc
%     Nseq = 'gggcgccc'; %ggg | ccc

function [LeftSide, RightSide] = divideNregion(Nseq,varargin)

%Convert Purines AG into X, and Pyrimidine CT into X.
clear Nclass;
Nclass(1:length(Nseq)) = ' ';
if length(varargin) == 2
    LeftNTs = varargin{1};
    RightNTs = varargin{2};
    
    %Need to insert the "|"
    LeftNT(1:length(LeftNTs)*2-1) = '|';
    LeftNT(1:2:end) = LeftNTs;
    RightNT(1:length(RightNTs)*2-1) = '|';
    RightNT(1:2:end) = RightNTs;
    
else
    LeftNT = 'a|g';
    RightNT = 'c|t';
end

PurLoc = regexpi(Nseq,LeftNT);
PyrLoc = regexpi(Nseq,RightNT);
% PurLoc = regexpi(Nseq,'c|t');
% PyrLoc = regexpi(Nseq,'a|g');
Nclass(PurLoc) = 'X';
Nclass(PyrLoc) = 'Y';

%Determine the cut point
CutMat1 = zeros(size(Nclass));
CutMat1(PurLoc) = 1;
CutMat1(PyrLoc) = -1;
SumCutMat1 = cumsum(CutMat1);

CutMat2 = zeros(size(Nclass));
CutMat2(PurLoc) = -1;
CutMat2(PyrLoc) = 1;
SumCutMat2 = fliplr(cumsum(fliplr(CutMat2)));

Max1 = max(SumCutMat1);
Max2 = max(SumCutMat2);

if Max1 <= 1 && Max2 <= 1
    LeftSide = '';
    RightSide = '';
    
elseif Max1 >= Max2
    CutLoc = find(Max1 == SumCutMat1);
    %If there is a tie, break by looking at local composition between the 2
    %peaks
    if length(CutLoc) > 1
        CutLocs = CutLoc(end-1:end);
        RightScore = SumCutMat2(CutLocs(1)+1:CutLocs(2)-1);        
        CutLoc2 = find(RightScore == max(RightScore));
        CutLoc = CutLocs(1) + CutLoc2(1)-1;
        
%         CutLoc = round(mean(CutLoc));
%         if mod(length(Nseq),2) == 1
%             if isempty(regexpi(Nseq(CutLoc),LeftNT))
%                 CutLoc = CutLoc + 1; %
%             end
%         end
    end            

    
    LeftSide = Nseq(1:CutLoc);
    RightSide = Nseq(CutLoc+1:end);    
    
elseif Max2 > Max1
    CutLoc = find(Max2 == SumCutMat2)-1; %Beecause we will cut to the "left side", need a -1.
    
    %If there is a tie, break by average
    if length(CutLoc) > 1
        CutLocs = CutLoc(1:2);
        LeftScore = SumCutMat1(CutLocs(1)+1:CutLocs(2));        
        CutLoc2 = find(LeftScore == max(LeftScore));
        CutLoc = CutLocs(1) + CutLoc2(end);
%         CutLoc = round(mean(CutLoc));
%         if mod(length(Nseq),2) == 1
%             if isempty(regexpi(Nseq(CutLoc),LeftNT))
%                 CutLoc = CutLoc + 1; %
%             end
%         end
    end    
    
    LeftSide = Nseq(1:CutLoc);
    RightSide = Nseq(CutLoc+1:end);    
end



if isempty(LeftSide)
    LeftSide = '';
elseif isempty(RightSide)
    RightSide = '';
end
