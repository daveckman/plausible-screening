function print_screening_results(alg_string, discrep_string, subset_indicators)
% Print screening results to screen

fprintf('Screening with %s ', alg_string)
if strcmp(alg_string,'PO') == 1 || strcmp(alg_string,'PO relaxed')
    fprintf('(%s)', discrep_string)
end
fprintf(': # solutions in subset = %d.\n', sum(subset_indicators))

end