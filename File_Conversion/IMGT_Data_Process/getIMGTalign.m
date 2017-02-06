function [Vseq, Vref, Vlen, Del5, Del3, Vname] = getIMGTalign(FID,SampName,Flocation)
%Go to the file location
fseek(FID,Flocation,'bof');

SampName2 = strrep(SampName,'_',' '); %IMGT is inconsistent and can look for space instead of underscore

TopSeq = '';
BotSeq = '';
BotSeqName = '';
CutLoc = 0;
while 1
    TextLine= fgetl(FID);
    if sum(regexp(TextLine,'\d+\.\s*','once')==1) > 0; break; end
    Check1 = ~isempty(strfind(TextLine,SampName) > 0);
    Check2 = ~isempty(strfind(TextLine,SampName2) > 0);
    if  Check1 || Check2 
        %Obtain the sample information
        if CutLoc == 0
            if Check1
                TextLine = strrep(TextLine,SampName,repmat(' ',1,length(SampName)));
            elseif Check2
                TextLine = strrep(TextLine,SampName2,repmat(' ',1,length(SampName2)));
            end
            CutLoc = regexp(TextLine,'[^\s]','once');
        end
        TopSeq = [TopSeq TextLine(CutLoc:end)];                        

        %Obtain the germline information
        TextLine = fgetl(FID);
        if isempty(BotSeqName)
            BotSeqName = TextLine(1:CutLoc-1);
        end
        BotSeq = [BotSeq TextLine(CutLoc:end)];
    end
end

%Determine the IMGT name for the gene
NameSeg = regexp(BotSeqName,'\s','split');
for j = 1:length(NameSeg)
    if sum(regexpi(NameSeg{j},'IGH[VDJ]')) > 0
        Vname = NameSeg{j};
        break
    end
end
ClassType = regexpi(Vname,'IGH([VDJ])','tokens');
ClassType = ClassType{1}{1};

%Determine start and end of sequences
switch upper(ClassType)
    case 'V'
        StartLoc = regexp(TopSeq,'\w','once');
        EndLoc = regexp(BotSeq,'\-');
        EndLoc = EndLoc(end);
    case 'D'
        StartEndLoc = regexp(BotSeq,'\-');
        StartLoc = StartEndLoc(1);
        EndLoc = StartEndLoc(end);
    case 'J'
        StartLoc = regexp(BotSeq,'\-','once');
        EndLoc = regexp(TopSeq,'\w');
        EndLoc = EndLoc(end);
end

%Remove the triple '...' located in the iMGT sequences.
DotLocs = BotSeq == '.';
DotLocs(1:StartLoc) = 0;
DotLocs(EndLoc:end) = 0;
TopSeq(DotLocs) = [];
BotSeq(DotLocs) = [];
EndLoc = EndLoc - sum(DotLocs);

%Find the deleted nt in the top seq, which must be added
TopDotLocs = find(TopSeq == '.');
TopSeq(TopDotLocs) = BotSeq(TopDotLocs);

%Remove the dots in the bot seq
CharLoc = regexp(BotSeq,'\w');
LeftCharLoc = CharLoc(CharLoc < StartLoc);
RightCharLoc = CharLoc(CharLoc > EndLoc);
if isempty(LeftCharLoc)
    Del5 = 0;
else
    Del5 = LeftCharLoc(end) - LeftCharLoc(1) + 1;
end
if isempty(RightCharLoc)
    Del3 = 0;
else
    Del3 = RightCharLoc(end) - RightCharLoc(1) + 1;
end
Vlen = EndLoc - StartLoc + 1;

Vseq = TopSeq(StartLoc:EndLoc); %(StartLoc:end);
Vref = BotSeq(StartLoc:EndLoc); %(StartLoc:end);



