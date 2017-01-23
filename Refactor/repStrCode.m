function S = repStrCode(S,VAR1,VAR2)
PatLoc = regexp(S,VAR1);
for j = 1:length(PatLoc)
    %Performing left side check
    ValidLeft = 0;
    if PatLoc(j) == 1; 
        ValidLeft = 1;
    elseif S(PatLoc(j)-1) == '.' %Exception - can't have a ".". Takes precedence over next clause.
        ValidLeft = 0;
    elseif ~isempty(regexp(S(PatLoc(j)-1),'[\W\s]+','once'))
        ValidLeft = 1;
    end
    
    %Performing right side check
    if ValidLeft == 1
        ValidRight = 0;
        if PatLoc(j)+length(VAR1)-1 == length(S);
            ValidRight = 1;
        elseif S(PatLoc(j)+length(VAR1)) == '.' %Exception - can't have a ".". Takes precedence of next clause.
            ValidRight = 0;
        elseif ~isempty(regexp(S(PatLoc(j)+length(VAR1)),'[\W\s]+','once'));
            ValidRight = 1;
        end
    end
    
    if ValidRight * ValidLeft == 0
        PatLoc(j) = 0;
    end
end
PatLoc(PatLoc == 0) = [];

%Iterative replacement
for j = 1:length(PatLoc)
    S1 = S(1:PatLoc(j)-1);
    S2 = S(PatLoc(j)+length(VAR1):end);
    S = [S1 VAR2 S2];
    PatLoc(j+1:end) = PatLoc(j+1:end) + length(VAR2) - length(VAR1);
end
