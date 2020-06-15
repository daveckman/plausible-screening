% (s,S) inventory problem
oracle_string = 'sS_oracle';
oracle_n_rngs = 1;
[A,B] = meshgrid(10:80,10:100);
feas_region = [A(:),B(:)];
feas_region = feas_region(feas_region(:,1)+14.5 < feas_region(:,2),:);
scrXn = [feas_region;...
    repmat(feas_region(feas_region(:,1)-feas_region(:,2)==-15,:),[5,1]);
    repmat(feas_region(feas_region(:,1)==10,:),[5,1]);
    repmat(feas_region(feas_region(:,1)==80,:),[5,1]);
    repmat(feas_region(feas_region(:,2)==10,:),[5,1]);
    repmat(feas_region(feas_region(:,2)==100,:),[5,1])];
k = 25;

% Reproduce same exp set
kmeans_rng = RandStream.create('mlfg6331_64');
opts = statset('Streams',kmeans_rng,'UseSubstreams',1);
[IDX, C] = kmeans(scrXn,25,'Options',opts);
exp_set = round(C);
exp_set(12,:) = [55,76]; % avoid singular covariance matrix my perturbing solution

n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'CRN'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant
clear('A', 'B', 'C', 'IDX', 'scrXn');