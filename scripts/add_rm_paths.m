function add_rm_paths(add_rm_flag)
% Add or remove necessary folders to the MATLAB path

folders = strcat({'..\'}, {'data','plots','problems','scripts','src','tools'});

switch add_rm_flag
    case 'add'
        for k = 1:length(folders)
            addpath(folders{k})
        end

    case 'remove'
        for k = 1:length(folders)
            rmpath(folders{k})
        end

    otherwise
        disp('\nERROR: Invalid argument for add_rm_paths(). Argument must be "add" or "remove".\n')
end

