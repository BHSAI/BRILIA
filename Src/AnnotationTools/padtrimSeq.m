%padtrimSeq will either pad or trim sequences so that SeqA and SeqB are the
%same length, and centered about an anchoring letter. Uses 'X' for padding.
%
%  NewSeqA = padtrimSeq(SeqA,AnchorA,LeftOpt,RightOpt)
%
%  [NewSeqA, NewSeqB] = padtrimSeq(SeqA,SeqB,AnchorA,AnchorB,LeftOpt,RightOpt)
%
%  [NewSeqA, NewSeqB, Aadj, Badj] = padtrimSeq(SeqA,SeqB,AnchorA,AnchorB,LeftOpt,RightOpt)
%
%  INPUT
%    SeqA/B: Sequences that you want to match with each other
%    AnchorA/B: Integer position # where you want SeqA/B to be same
%    LeftOpt/RightOpt: Option to set left and right lengths
%      'Min'    set to same min length with respect to anchor point
%      'Max'    set to same max length with respect to anchor point
%      'A'      set to same length with repsect to SeqA and anchor point
%      'B'      set to same length with repsect to SeqB and anchor point
%      # int    set to same N number of letters left/right of anchor point
%
%  OUTPUT
%    SeqA/B: new sequences that have the same-length sequences after
%      padding and/or trimming
%    Aadj/Badj: 1x2 matrix showing how many letters were removed (- number)
%      or added (+ number)
%  
%  EXAMPLE
%    SeqA = 'ACGTGCGTA'
%    SeqB = 'CGAGTGCATGTG'
%    AnchorA = 4;
%    AnchorB = 5;
%
%    [NewSeqA, Aadj] = padtrimSeq(SeqA,AnchorA,2,7)
%    NewSeqA = 
%       CGTGCGTAXX
%    Aadj = 
%       -1    2
%
%    [NewSeqA, NewSeqB] = padtrimSeq(SeqA,SeqB,AnchorA,AnchorB,'min','max')
%    NewSeqA = 
%       ACGTGCGTAXX
%    NewSeqB = 
%       GAGTGCATGTG
%
%    [NewSeqA, NewSeqB, Aadj, Badj] = padtrimSeq(SeqA,SeqB,AnchorA,AnchorB,2,9)
%    NewSeqA = 
%       CGTGCGTAXXXX
%    NewSeqB = 
%       AGTGCATGTGXX
%    Aadj =
%       -1    4 
%    Badj =
%       -2    2
%
%  See also findVDJmatch, alignSeq, findGeneMatch, getGeneSeed

function varargout = padtrimSeq(varargin)
%Determine if you have a single or double sequence input
InvalidInput = 0;
if length(varargin) == 4 %Single input
    [SeqA, AnchorA, LeftOpt, RightOpt] = deal(varargin{:});   
    InputCount = 1;
    if isempty(SeqA);
        InvalidInput = 1;
    end
else %Two seuqence input
    [SeqA, SeqB, AnchorA, AnchorB, LeftOpt, RightOpt] = deal(varargin{:});
    InputCount = 2;
    if isempty(SeqA) || isempty(SeqB)
        InvalidInput = 1;
    end
end

if InvalidInput == 1
    %Format optional outputs
    if nargout >= 1
        varargout{1} = SeqA;
        if nargout >= 2
            varargout{2} = SeqB;
            if nargout >= 3
                varargout{3} = [0 0];
                if nargout >= 4
                    varargout{4} = [0 0];
                end
            end
        end
    end
    return
end

if InputCount == 1
    %Setup default adjust
    Aadj = [0 0];

    %Determine starting stats
    LeftA = AnchorA - 1;
    RightA = length(SeqA) - AnchorA;

    %Determine final Left and Right side lengths
    SetLeft = LeftOpt;
    SetRight = RightOpt;

    %Determine what needs to be adjusted
    Aadj(1) = -LeftA + SetLeft;
    Aadj(2) = -RightA + SetRight;

    %Create the padded/trimmed SeqA and SeqB
    if Aadj(1) > 0
        LeftApad = repmat('X',1,Aadj(1));
        S1a = 1;
    else
        LeftApad = '';
        S1a = -Aadj(1) + 1;
    end

    if Aadj(2) > 0
        RightApad = repmat('X',1,Aadj(2));
        S2a = length(SeqA);
    else
        RightApad = '';
        S2a = length(SeqA) + Aadj(2);
    end

    NewSeqA = [LeftApad SeqA(S1a:S2a) RightApad];

    %Format optional outputs
    if nargout >= 1
        varargout{1} = NewSeqA;
        if nargout >= 2
            varargout{2} = Aadj;
        end
    end
    
elseif InputCount == 2
    %Setup default adjust
    Aadj = [0 0];
    Badj = [0 0];

    %Determine starting stats
    LeftA = AnchorA - 1;
    RightA = length(SeqA) - AnchorA;

    LeftB = AnchorB - 1;
    RightB = length(SeqB) - AnchorB;

    %Determine final Left side length
    if isnumeric(LeftOpt)
        SetLeft = LeftOpt;
    else
        switch lower(LeftOpt)
            case 'min'
                SetLeft = min([LeftA LeftB]);
            case 'max'
                SetLeft = max([LeftA LeftB]);
            case 'a'
                SetLeft = LeftA;
            case 'b'
                SetLeft = LeftB;
        end
    end

    %Determine final Right side length
    if isnumeric(RightOpt)
        SetRight = RightOpt;
    else
        switch lower(RightOpt)
            case 'min'
                SetRight = min([RightA RightB]);
            case 'max'
                SetRight = max([RightA RightB]);
            case 'a'
                SetRight = RightA;
            case 'b'
                SetRight = RightB;
        end
    end

    %Determine what needs to be adjusted
    Aadj(1) = -LeftA + SetLeft;
    Aadj(2) = -RightA + SetRight;
    Badj(1) = -LeftB + SetLeft;
    Badj(2) = -RightB + SetRight;

    %Create the padded/trimmed SeqA and SeqB
    if Aadj(1) > 0
        LeftApad = repmat('X',1,Aadj(1));
        S1a = 1;
    else
        LeftApad = '';
        S1a = -Aadj(1) + 1;
    end

    if Aadj(2) > 0
        RightApad = repmat('X',1,Aadj(2));
        S2a = length(SeqA);
    else
        RightApad = '';
        S2a = length(SeqA) + Aadj(2);
    end

    if Badj(1) > 0
        LeftBpad = repmat('X',1,Badj(1));
        S1b = 1;
    else
        LeftBpad = '';
        S1b = -Badj(1) + 1;
    end

    if Badj(2) > 0
        RightBpad = repmat('X',1,Badj(2));
        S2b = length(SeqB);
    else
        RightBpad = '';
        S2b = length(SeqB) + Badj(2);
    end

    NewSeqA = [LeftApad SeqA(S1a:S2a) RightApad];
    NewSeqB = [LeftBpad SeqB(S1b:S2b) RightBpad];

    %Format optional outputs
    if nargout >= 1
        varargout{1} = NewSeqA;
        if nargout >= 2
            varargout{2} = NewSeqB;
            if nargout >= 3
                varargout{3} = Aadj;
                if nargout >= 4
                    varargout{4} = Badj;
                end
            end
        end
    end
end

