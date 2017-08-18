%formatSeq takes a Seq and Classifier, and adds spaces between each class,
%lower case mismatches, upper case matches.
%
%  FormattedSeq = formatSeq(Seq,Classifier)
%
%  EXAMPLE 
%    Seq = 'ACGTTGTGCACGTGTTTATAT';
%    Classifier = 'VVVvVVMMMDdDDNNJJjJJj';
%    FormattedSeq = formatSeq(Seq,Classifier);
%    FormattedSeq = 
%           ACGtTG tgc AcGT gt TTaTAt

function FormattedSeq = formatSeq(Seq,Classifier)
%Make sure no input is a cell
if iscell(Seq)
    Seq = Seq{1};
end
if iscell(Classifier)
    Classifier = Classifier{1};
end

%Look for the lower cases and format those NTs to lower case
MisMatchLoc = isstrprop(Classifier,'lower');
Seq = upper(Seq);
Seq(MisMatchLoc) = lower(Seq(MisMatchLoc));

%Convert VD and DJ p-nucleotides into M or N respecitivley, and make NTs lower case.
Classifier(regexpi(Classifier,'b')) = 'm';
Classifier(regexpi(Classifier,'p')) = 'n';
MmatchLoc = regexpi(Classifier,'m');
NmatchLoc = regexpi(Classifier,'n');
Seq([MmatchLoc NmatchLoc]) = lower(Seq([MmatchLoc NmatchLoc]));

ClassUnique = 'VMDNJ';

%Assemble the formatted sequence.
FormattedSeq(1:length(Seq)+length(ClassUnique)-1) = ' ';
q = 1; 
for j = 1:length(ClassUnique)
    ClassLoc = regexpi(Classifier,ClassUnique(j));
    FormattedSeq(q:q+length(ClassLoc)-1) = Seq(ClassLoc);
    q = q+length(ClassLoc)+1;
end

