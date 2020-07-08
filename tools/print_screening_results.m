function print_screening_results(problem_string, feas_region, exp_set, alg_string, discrep_string, fn_props, subset_indicators)
% Print screening results to screen

fprintf('\nResults for problem "%s"\n', problem_string)
fprintf('\t# of solutions in feasible region: \t\t%d\n', size(feas_region,1))
fprintf('\t# of solutions in experimental set: \t%d\n', size(exp_set,1))
fprintf('__________________________________________________________\n\n')
fprintf('Screening with %s\n', alg_string)
if strcmp(alg_string,'PO') == 1 || strcmp(alg_string,'PO relaxed')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\tfunctional property: \t\t\t\t\t%s\n', fn_props)
end
fprintf('\t# of solutions in subset: \t\t\t\t%d\n', sum(subset_indicators))

end