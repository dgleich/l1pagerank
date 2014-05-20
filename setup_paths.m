if ispc
    
    mypath=fileparts(mfilename('fullpath'));
    addpath(mypath)
    addpath('e:\dev\matlab-bgl\4.0');
    addpath('e:\doc\My Dropbox\matlab-extern\cvx\');
    cvx_setup
    
else
    
    mypath=fileparts(mfilename('fullpath'));
    addpath(mypath);
    addpath('~/dev/mcode');
    addpath('~/dev/matlab-bgl');
    addpath('~/Dropbox/matlab-extern/cvx/');
    cvx_setup
    
    if ismac
        
        cwd = pwd;
        try
            cd /Library/gurobi560/mac64/matlab/
            gurobi_setup
        catch
            cd(cwd)
        end
        cd(cwd)
        setenv('PYTHON','python2.7');
    else
        fprintf('Did you use ~/.gurobi_setup????');
        addpath('/opt/gurobi560/linux64/matlab/');
        gurobi_setup;
        
    end

end

addpath(fullfile(mypath,'metis-5.0.2','metismex'));