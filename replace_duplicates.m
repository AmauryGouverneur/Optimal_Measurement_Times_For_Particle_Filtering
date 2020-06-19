function newLine = replace_duplicates(line,accessible_meas)
    newLine = line;
    duplicatesIndices = logical([0 ~diff(line)]);
    yetIn = line(~duplicatesIndices);
    addable = setdiff(accessible_meas(randperm(length(accessible_meas))),yetIn,'stable');
    newLine(duplicatesIndices) = addable(1:sum(duplicatesIndices));
end
