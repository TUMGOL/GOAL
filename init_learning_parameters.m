%(c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
%Muenchen, 2012. Contact: simon.hawe@tum.de
function para = init_learning_parameters()
    para.p        = 0;
    para.q        = 2;
    para.kappa    = 1e5;
    para.nu       = 1e4;
    para.mu    = 1e6;
    para.max_iter = 300;
    % LogAbs, LogSquare, PNormAbs, PNormSquare, AtanAbs, AtanSquare
    para.Sp_type  = 'LogSquare'; 
    para.Omega    = [];
    para.verbose  = 1;
    para.logger   = [];
end