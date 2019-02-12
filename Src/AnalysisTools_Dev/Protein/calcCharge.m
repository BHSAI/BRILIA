%calcCharge will calculate the net charge of an amino acid sequence at a
%certain pH, and with or without the C- and N- terminus. The pKa values
%used are from EMBOSS package.
%
%  Charge = calcCharge(Seq)
%
%  Charge = calcCharge(Seq, Param, Value, ...)
%
%  Info = calcCharge()
%
%  INPUT
%    Seq: string or cell of string of amino acid sequences
%    *If empty, will return the Info structure
%
%    Param        Value    Details
%    ------------ -------- --------------------------------------
%    pH           #        pH to evaluate the protein charge (Default 7.2)
%    IncludeTerm  1 or 0   1 = include C and N terminus (Default)
%                          0 = exclude C and N terminus
%    Info         Nx3 cell New values of pKa to use. 
%                          Col1 = AA letter and 'Nterm', 'Cterm'
%                          Col2 = pKa of each entity
%                          Col3 = charge +1 or -1 
%
%  OUTPUT
%    Charge: value or matrix of net charges
%    Info: structure of HPI data
% 
%  EXAMPLE
%    Seq = 'CARDYLWL';
%    NetCharge = calcCharge(Seq, 'IncludeTerminus', false)
%    NetCharge =
%        -0.0485
function NetCharge = calcCharge(Seq, varargin)
P = inputParser;
addParameter(P, 'pH', 7.2, @isnumeric);
addParameter(P, 'IncludeTerminus', true, @(x) islogical(x) || x == 0 || x == 1);
addParameter(P, 'Info', [], @(x) isempty(x) || (iscell(x) && size(x, 1) == 3));
parse(P, varargin{:})
IncludeTerminus = P.Results.IncludeTerminus;
pH = P.Results.pH;
Info = P.Results.Info;

if ~isempty(Info) %Use custom pKa AND charge of the conjugate acid
    Code = Info(:, 1);
    pKa  = cell2mat(Info(:, 2));
    Charge = cell2mat(Info(:, 3));
else  %Store the EMBOSS pKa AND charge of the conjugate acid
    Code   = {'Nterm' 'Cterm'    'C'    'D'    'E'    'H'     'K'     'R'     'Y'}';
    pKa    = [   8.60    3.60   8.50   3.90   4.10   6.50   10.80   12.50   10.10]';
    Charge = [      1      -1     -1     -1     -1      1       1       1      -1]';
end

if nargin == 0
    InfoCell = [Code num2cell(pKa)]';
    NetCharge = struct(InfoCell{:});
    return
end

Seq = upper(Seq);

Count = zeros(size(Code));
if IncludeTerminus
    Count(1:2) = 1;
end

if iscell(Seq)
    NetCharge = zeros(length(Seq), 1);
    for k = 1:length(Seq)
        for j = 3:length(Code)
            Count(j) = sum(Seq{k} == Code{j});
        end
        NetCharge(k) = calcNetCharge(pH, Count, pKa, Charge);
    end
else
    for j = 3:length(Code)
        Count(j) = sum(Seq == Code{j});
    end
    NetCharge = calcNetCharge(pH, Count, pKa, Charge);    
end

function NetCharge = calcNetCharge(pH, Count, pKa, Charge)
PosLoc = Charge > 0;
NegLoc = Charge < 0;
T = 10.^(pH - pKa);
NetCharge =  sum(Count(PosLoc) .* 1./(1 + T(PosLoc))) - ...
             sum(Count(NegLoc) .* T(NegLoc)./(1 + T(NegLoc)));