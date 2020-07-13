function print_screening_results(alg_string, discrep_string, subset_indicators)
% Print screening results to screen

fprintf('\nScreening with %s\n', alg_string)
if strcmp(alg_string,'PO') == 1 || strcmp(alg_string,'PO relaxed')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
end
fprintf('\t# of solutions in subset: \t\t\t\t%d\n\n', sum(subset_indicators))

end