%% External Loading of 1D Lumped Parameter Muscle Models

clear all
addpath('..\');
addpath('..\..\..\1D-Hill-2-Component-Model\code');
MEMBUFF = 20000;

%% 
% Global plotting options
fnum = 1;
opt_grid = 'on';
opt_hold =  'on';
splotx = 2;
sploty = 2;

%% Global muscle properties
% $L_{rest}$ is the muscle resting length, $F_{max}$ is the maximum force the 
% muscle can produce, $F_{max}^{ECC}$ is the maximal force (parameter) under eccentrinc 
% loading and $V_{max}$ is the maximum velocity of shortening.

% Minimum length 
L_MIN = 0.2; %mm 
% Maximum length 
L_MAX = 0.8; %mm
L_REST = 0.5; %mm
F_MAX = 5; %N
V_MAX = -1.5; %mm/s
F_MAXECC = 2*F_MAX; %N

%% 
% Consider an external loading as follows:
% 
% # Stretch the muscle passively to some intiial length: $L_{INIT}$
% # Activate the muscle to some activity level: $\alpha_{FIXED}$ (for now, consider 
% full activation) to give an initial (isometric) force: $F_{INIT}$
% # Free one end of the muscle and instantly apply an external load: $F_{EXT}$
% # Find the resulting final length of the muscle

L_INIT = 0.625;
alpha = 1;
V = 0; % for an isometric contraction

% Find total muscle force
[F_MUSC_INIT,F_ACT_INIT,F_PASS_INIT] = force_muscle(L_INIT, ... 
        L_REST, V_MAX, V, F_MAX, alpha, F_MAXECC);   
    
%% 
% Consider that the muscle has to undergo a dynamic contraction from $L_{INIT}$ 
% to $L_{FIN}$ (while contracting against $F_{EXT}$). The velocity of this contraction 
% is given by the force-velocity law; then the first step is to see if the external 
% load would cause an eccentric or concentric contraction. 
% 
% Find the velocity for the given external loading (assume for now that it 
% is below the maximal eccentric force. Given that any number of force-velocity 
% laws of the form $F_{VEL} = f(V,V_{MAX}, ...)$ can be used, and there may not 
% exists a analytical solution in the form $V = f(F,V_{MAX}, ...)$. Therefore 
% solve iteratively for a given external force to find the corresponding velocity.
% 
% Once the velocity of contraction is found _at the current muscle length_, 
% the next steps are:
% 
% # Update the muscle length
% # Update the force-velocity relationship
% # Find the new contraction velocity (for the fixed $F_{EXT}$)
% 
% What is expected is that the velocity of contraction tends towards zero 
% - and the length is found where $F_{EXT}$ is the  $F_{MAX}^{LOCAL}$ - and the 
% muscle is at rest.
%
% *Remark*: Here the acceleration and deceleration phases of the muscle are
% not explicitly modeled. In that the muscle has no mass. Any acceleration
% or deceleration that arises is due to the dynamic change in the
% force-velocity relationship; i.e. $F_{VEL} = f(V,L,V_{MAX},F_{MAX}^{LOCAL})$
%
% 
% *Force enhancement* occurs for the case when a fully activated muscle
% is undergoing an eccentrinc contraction; infact only on the decending
% limb of the force-length response. Therefore, the computation steps
% are:
% 
% # Check if current muscle length L_{TOT} is on the decending limb by
% taking a small perturbation in the length. If on decending limb,
% continue, if not escape.
% # Check if the velocity will cause an elongation of the muscle
% # Determine the change in length; based on the velocity & time-step
% # Use the force-enhancement law to compute the resulting force
% enhancement F_{ENH}
% 
% This is then added to the total muscle force such that:
%
% $$F_{MUSC} = F_{ACT} + F_{PASS} + F_{ENH}$$
%

% Specify the external load
F_EXT = 0.5; 
% define a time-step
tStep = 0.01; %s
tol = 0.05;

V_C = ones(1,MEMBUFF);
L_CURR(MEMBUFF) = zeros;
F_MUSC(MEMBUFF) = zeros;
F_PASS(MEMBUFF) = zeros;
F_ACT(MEMBUFF) = zeros;
F_ENH(MEMBUFF) = zeros;
i = 1;
L_CURR(i) = L_INIT;
V_C(i) = 0;
V_PREV = 1;
L_PREV = L_MAX/2;
F_DIFF = 1;
F_PASS(i) = F_PASS_INIT;
ENH_START = 0;

RFH = 1;

while abs(F_DIFF) > tol
%while abs(V_PREV) > tol && L_PREV > L_MIN && L_PREV < L_MAX
    % find $F_{MAX}^{LOCAL}$
    % find the "local" F_{MAX}^{LOCAL}, i.e. the F_{MAX} at the given (inital) length L_{INIT}
    F_MAXLOCAL = force_length(L_CURR(i), L_REST);
    F_MAXLOCAL = F_MAXLOCAL*F_MAX;

    
    % eccentric on dec. & asc. limb; then residual force enhancement
    if RFH == 1
        if V_C(i) > 0
            F_MAXLOCAL = F_MUSC(i-1);
        end
    end
    
    % contractile velocity for a given F_{EXT}
    % remove the passive component
    F_EXT_ACC = F_EXT - F_PASS(i);
    V_C(i) = contractile_vel(F_EXT_ACC, V_MAX, F_MAXECC, F_MAXLOCAL, L_CURR(i), L_REST);
    % update muscle length
    L_CURR(i+1) = L_CURR(i) + V_C(i)*tStep;
    %V_PREV = V_C(i);
    %L_PREV = L_CURR(i); 
    
 
    % Find the total muscle force (without RESIDUAL FORCE ENHACENEMENT)
    % at zero velocity i.e. here F_{MAX}^{LOCAL} = F_{EXT}
    [F_MUSC(i),F_ACT(i),F_PASS(i)] = force_muscle(L_CURR(i), ... 
        L_REST, V_MAX, 0, F_MAX, alpha, F_MAXECC); 
    
    % RESIDUAL FORCE ENHANCEMENT
    if RFH == 1
        if V_C(i) < 0
            % do nothing
        elseif V_C(i) > 0 && L_CURR(i) > L_REST
            if ENH_START == 0
                L_ENH_START = L_CURR(i);
                F_ENH_START = F_MUSC(i);
            end
            ENH_START = 1;
            % residual force enhancement at current length
            F_ENH(i) = force_enhance(L_CURR(i), F_MUSC(i), L_ENH_START, F_ENH_START);
        end
    end
    
    % ADD IN RESIDUAL FORCE ENHANCEMENT
    %F_MUSC(i) = F_MUSC(i)+F_ENH(i);
    F_MUSC(i) = max(F_MUSC(i),F_MUSC(i)+F_ENH(i));
    
    % for convergence
    F_DIFF = F_MUSC(i) - F_EXT;

    % update increment
    i = i+1;      

end


% time vector (for plotting)
tvec = 0:tStep:tStep*(i-2);

xvec = tvec;
yvec = L_CURR(1:i-1);
ftitle = 'Length of Muscle';
xtitle = 't (s)';
ytitle = 'Length (mm)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)


fnum = fnum + 1;
xvec = tvec;
yvec = V_C(1:i-1);
ftitle = 'Velocity of Muscle';
xtitle = 't (s)';
ytitle = 'Velocity (mm/s)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)


fnum = fnum + 1;

xvec = tvec;
yvec = F_MUSC(1:i-1);
ftitle = 'Total Muscle Force';
xtitle = 't (s)';
ytitle = 'F_{MUSC} (N)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

fnum = fnum + 1;
xvec = L_CURR(1:i-1);
yvec = F_MUSC(1:i-1);
ftitle = 'Component Muscle Forces';
xtitle = 'L (mm)';
ytitle = 'F_{MUSC} (N)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

