%dispList will display all contents of a string cell, listing it by number
%too. Used mainly for interfacing with the user what option to select.
%  
%  dispList(A)
%
%  EXAMPLE 
%    A = {'What'; 'Is'; 'This'; 1322 ; 2332; 31423}
%    dispList(A)
%    A = 
%         1) What
%         2) Is
%         3) This
%         4) 1322
%         5) 2332
%         6) 31423

function dispList(A)
if ischar(A)
    A = {A};
end

if min(size(A)) ~= 1
    error('Needs to be a 1xN or Nx1 cell matrix')
end

for j = 1:length(A)
    if ischar(A{j})
        disp(sprintf('%d) %s',j,A{j}));
    elseif isnumeric(A{j})
        disp(sprintf('%d) %s',j,num2str(A{j})));
    end
end
