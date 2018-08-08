function F_ENH = force_enhance(L_TOT, F_MUSC, L_ENH_START, F_ENH_START)
    % there should be no enhancement on the f-l curve
    gamma = 10;
    F_ENH_TOT = gamma*(L_TOT-L_ENH_START)+F_ENH_START;
    
    F_ENH = F_ENH_TOT - F_MUSC;


