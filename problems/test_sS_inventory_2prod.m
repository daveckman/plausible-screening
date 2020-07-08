% (s,S) inventory problem with 2 products
oracle_string = 'sS_2prod_oracle';
oracle_n_rngs = 1;
[A,B] = meshgrid(10:80,10:100);
single_feas_region = [A(:),B(:)];
single_feas_region = single_feas_region(single_feas_region(:,1)+14.5 < single_feas_region(:,2),:);

% Cross designs
[a, b] = ndgrid(1:size(single_feas_region,1), 1:size(single_feas_region,1));
feas_region = [single_feas_region(a,:), single_feas_region(b,:)];

scrXn = [single_feas_region;...
    repmat(single_feas_region(single_feas_region(:,1)-single_feas_region(:,2)==-15,:),[5,1]);
    repmat(single_feas_region(single_feas_region(:,1)==10,:),[5,1]);
    repmat(single_feas_region(single_feas_region(:,1)==80,:),[5,1]);
    repmat(single_feas_region(single_feas_region(:,2)==10,:),[5,1]);
    repmat(single_feas_region(single_feas_region(:,2)==100,:),[5,1])];
k = 25;

% Reproduce same exp set
kmeans_rng = RandStream.create('mlfg6331_64');
opts = statset('Streams',kmeans_rng,'UseSubstreams',1);
[IDX, C] = kmeans(scrXn,25,'Options',opts);
single_exp_set = round(C);
single_exp_set(12,:) = [55,76]; % avoid singular covariance matrix my perturbing solution

% Cross designs
[a, b] = ndgrid(1:size(single_exp_set,1), 1:size(single_exp_set));
exp_set = [single_exp_set(a,:), single_exp_set(b,:)];

n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}
clear('A', 'B', 'C', 'IDX', 'scrXn', 'single_feas_region', 'single_exp_set', 'a', 'b');