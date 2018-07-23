%dispList will display all contents of a string cell, listing it by number
%too. Used mainly for interfacing with the user what option to select.
%  
%  dispList(List)
%
%  INPUT
%    List: Mx1 cell list of items to show
%
%  EXAMPLE 
%    List = {'What'; 'Is'; 'This'; 1322 ; 2332; 31423}
%    dispList(List)
%         1) What
%         2) Is
%         3) This
%         4) 1322
%         5) 2332
%         6) 31423

function dispList(List)
if ischar(List)
    List = {List};
end

if min(size(List)) ~= 1
    error('Needs to be a 1xN or Nx1 cell matrix')
end

for j = 1:length(List)
    if ischar(List{j})
        fprintf('%d) %s\n',j,List{j});
    elseif isnumeric(List{j})
        fprintf('%d) %s\n',j,num2str(List{j}));
    end
end
