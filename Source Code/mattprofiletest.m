clear all; close all; clc
m = 50;
test_sS_inventory;
[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);


profile on
[~, ~, ~, ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params);
profile off
profile viewer