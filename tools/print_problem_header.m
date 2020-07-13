function print_problem_header(problem_string, feas_region, exp_set, fn_props)
% Print problem header to screen

fprintf('\nResults for problem "%s"\n', problem_string)
fprintf('\t# of solutions in feasible region: \t\t%d\n', size(feas_region,1))
fprintf('\t# of solutions in experimental set: \t%d\n', size(exp_set,1))
fprintf('\tfunctional property: \t\t\t\t\t%s\n', fn_props)
fprintf('__________________________________________________________\n\n')

end