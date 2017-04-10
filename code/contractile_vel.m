function V_C = contractile_vel(F_EXT, V_MAX, F_MAXECC, F_MAXLOCAL)

% 1. Find if the external load causes ecc., con. or iso. contraction
F_RATIO = F_EXT/F_MAXLOCAL;

if F_RATIO < 1 % CONCENTRIC
    C_TYPE = -1;
elseif F_RATIO > 1 % ECCENTRIC
    C_TYPE = 1;
elseif F_RATIO == 1 % ISOMETRIC
    C_TYPE = 0;
end 

% find the corresponding velocity of contraction at given muscle length
F_VEL = -1;
tol = 0.001;
vStep = tol/10;
V = 0;
while abs(F_VEL-F_EXT/F_MAXLOCAL) > tol && abs(V) < abs(V_MAX)
    if C_TYPE < 0 % CONCENTRIC
        V = V-vStep;
    elseif C_TYPE > 0 % ECCENTRIC - DYNAMIC FORCE ENHANCEMENT
        V = V+vStep;
    elseif C_TYPE == 0 % ISOMETRIC
        V = 0;
    end
    F_VEL = force_vel(V_MAX, V, F_MAXLOCAL,F_MAXECC);
end
V_C = V; % velocity of contraction