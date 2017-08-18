%DUPLICATE of parseInput to remove depencies.
%
%Parse input will parse a cell, structure, or command-line input and return
%the results in either a structured input or cell input. Use this parsing
%for dealing with codes that can be either used in a command-line or Matlab
%environment, as it will take care:
%  1) leading dashes in front of parameter names [EX: '-saveas' => 'saveas]
%  2) lower/upper case mismatches [EX: 'SaveAS' => 'SaveAs']
%  3) evaluating numeric strings [EX: '1,2' => [1 2], '[1:3]' => [1 2 3].
%
%  [Ps, Pu, ReturnThis] = parseInputT(P, Varargin{:})
%
%  [Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInputT(P, Varargin{:})
%
%  INPUT
%    P: inputParser object
%    Varargin: 1xM varargin cell input to be parsed
%
%  OUTPUT
%    Ps: the used structured input P.(Param) = Value.
%    Pu: the unused structured inputs P.(UnusedParam) = Value.
%    ReturnThis: 1 if the function should just return stop due to
%      'getinput' or 'getvarargin' command, 0 if it should continue with
%      the rest of the code.
%    ExpPs: expanded version of Ps as a cell
%    ExpPu: expanded version of Pu as a cell
%
%  EXAMPLE USE IN CODE
%    function out = func(vararin);
%    P = inputParser
%    P = addParameter(P, 'Var', 1, @isnumeric);
%    [Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInputT(P, varargin{:})
%    if ReturnThis
%        varargout = {Ps, Pu, ExpPs, ExpPu};
%        return;
%    end
%    Var = P.Var; %This function can used the right variables
%    Out2 = funct2(ExpPu{:}); %Subfunctions might use the other inputs.
%    Out3 = funct3(Pu); %Subfunctions could use Pu directly

function [Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInputT(P, varargin)
P.CaseSensitive = 0;
P.KeepUnmatched = 1;
P.PartialMatching = 0;
P.StructExpand = 1;

varargin = cleanCommandLineInput(varargin{:});

Pu = [];
ExpPs = [];
ExpPu = [];
ReturnThis = 0;
if ~isempty(varargin) && ischar(varargin{1})
    if ismember(lower(varargin{1}), {'getinput', 'getinputs'})
        ReturnThis = 1;
        parse(P);
        Ps = P.Results;
        return;
    elseif ismember(lower(varargin{1}), {'getvarargin', 'extractvarargin'})
        ReturnThis = 1;
        varargin(1) = [];
    end
end
parse(P, varargin{:});
Ps = P.Results;
Pu = P.Unmatched;
ExpPs = expandStruct(Ps);
ExpPu = expandStruct(Pu);

%Converts a structured input to a cell input format
function CellInput = expandStruct(StructInput)
Fields = fieldnames(StructInput);
Values = struct2cell(StructInput);
CellInput = [Fields'; Values'];
CellInput = CellInput(:)';

%Tries to evalues numbers, and remove leading dashes in field names (eg,
%EX:  '-saveas' => 'saveas'
%EX:  '[1,3]'   =>  [1 3]
function varargin = cleanCommandLineInput(varargin)
%Basic trial-error interpretation of number inputs
for j = 1:length(varargin)
    if ischar(varargin{j}) && ~isempty(varargin{j})
        InputStr = varargin{j};
        if isempty(regexpi(InputStr, '[^0-9\,\.\:\[\]\(\)]', 'once'))
            BrackLoc = regexpi(InputStr, '\[\]\(\)');
            if ~(InputStr(1) == '.') %In case users inputs a '.3' prefix as a string, skip.
                InputStr(BrackLoc) = [];
                try
                    varargin{j} = eval(['[' InputStr ']']);
                    continue;
                catch
                end
            end
        end
        
        if InputStr(1) == '-'
            NonDashLoc = 1;
            while InputStr(NonDashLoc) == '-'
                NonDashLoc = NonDashLoc + 1;
            end
            varargin{j} = InputStr(NonDashLoc:end);
            continue;
        end
    end
end
