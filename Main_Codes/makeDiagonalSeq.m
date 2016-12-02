%makeDiagonalSeq will take a SeqA, and then return a character array with
%SeqA moving from right to left. This is good for  using convolveSeq.
%
%  EX:
%    Seq1 = 'ABC'
%    Seq2 = 'ABCDEF';
%    makeDiagonalSeq(Seq1,Seq2)     'ABCDEFGH','ABC')
%      ans =     
%         ---A
%         --AB
%         -ABC
%         ABCD
%         BCD-
%         CD--
%         D---
function A = makeDiagonalSeq(Seq1,Seq2,varargin)
Tlen = length(Seq1)+length(Seq2)-1;
SeqT = [repmat('-',Tlen,1);Seq2'];
RectSeq = repmat(SeqT,1,Tlen+1);
RectSeq = RectSeq(:);
RectSeq(end+1:end+Tlen+1) = '-';
A = reshape(RectSeq,length(RectSeq)/(Tlen+1),Tlen+1)';

if length(varargin) == 2
    s1 = varargin{1};
    s2 = varargin{2};
    A = A(s1:s2,end-length(Seq2)-length(Seq1)+1:end-length(Seq2));
else
    A = A(:,end-length(Seq2)-length(Seq1)+1:end-length(Seq2));
end    

% Tlen = size(RectSeq,1)+1;
% RectSeq = RectSeq(:);
% Tadd = Tlen - mod(length(RectSeq),Tlen);
% RectSeq(end+1:end+Tadd) = RectSeq(1:Tadd);
% Output = reshape(RectSeq,Tlen,length(RectSeq)/Tlen)';
% To = length(Seq1)-length(Seq2)+1;
% Output = Output(To+s1-1:To+s2-1,1:length(Seq2));
% % 
% % 
% % SeqT = [repmat('-',length(Seq1)-1,1); Seq1'];
% RectSeq = repmat(SeqT,1,length(SeqT));
% Tlen = size(RectSeq,1)+1;
% RectSeq = RectSeq(:);
% Tadd = Tlen - mod(length(RectSeq),Tlen);
% RectSeq(end+1:end+Tadd) = RectSeq(1:Tadd);
% Output = reshape(RectSeq,Tlen,length(RectSeq)/Tlen)';
% To = length(Seq1)-length(Seq2)+1;
% Output = Output(To+s1-1:To+s2-1,1:length(Seq2));
