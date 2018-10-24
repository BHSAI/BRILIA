%getAllHeaderVar combines getLightHeaderVar and getHeavyHeaderVar into a
%single function, and also returns the Chain ('H','L','HL') depending on
%the VDJheader that is used. This is used to simplify codes.
%
%  [H, L] = getAllHeaderVar(VDJheader)
%
%  [H, L, Chain] = getAllHeaderVar(VDJheader)
%
%  [H, L, Chain, N] = getAllHeaderVar(VDJheader)
%
%  [H, L, Chain, N, C] = getAllHeaderVar(VDJheader)
%
%  INPUT
%    VDJheader: 1xN cell of VDJdata header names
%
%  OUTPUT
%    H: structure containing heavy chain header locations
%    (getHeavyHeaderVar.m)
%    L: structure containing light chain header locations
%    (getLightHeaderVar.m)
%    Chain: 'H', 'L', or 'HL', depending on VDJheader
function [H, L, varargout] =  getAllHeaderVar(VDJheader)
warning('%s: This is deprecated. Use getVDJmapper.', mfilename);

if nargout <= 3
    H = getHeavyHeaderVar(VDJheader);
    L = getLightHeaderVar(VDJheader);
    if nargout >= 3 %Return the chain info too
        if H.SeqLoc > 0 && L.SeqLoc > 0
            varargout{1} = 'HL';
        elseif H.SeqLoc > 0
            varargout{1} = 'H';
        elseif L.SeqLoc > 0
            varargout{1} = 'L';
        else
            varargout{1} = 'none';
        end
    end
else
    [H, Nh, Ch] = getHeavyHeaderVar(VDJheader);
    [L, Nl, Cl] = getLightHeaderVar(VDJheader);
    if nargout >= 3 %Return the chain info too
        if H.SeqLoc > 0 && L.SeqLoc > 0
            varargout{1} = 'HL';
        elseif H.SeqLoc > 0
            varargout{1} = 'H';
        elseif L.SeqLoc > 0
            varargout{1} = 'L';
        else
            varargout{1} = 'none';
        end
        if nargout >= 4 %Return N
            varargout{2} = unique([Nh; Nl]);
            if nargout >= 5 %Return C
                varargout{3} = unique([Ch; Cl]);
            end
        end
    end
end
