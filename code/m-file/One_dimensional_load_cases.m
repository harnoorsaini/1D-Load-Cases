%% External Loading of 1D Lumped Parameter Muscle Models

clear all
addpath('..\');
addpath('..\..\..\1D-Hill-2-Component-Model\code');
MEMBUFF = 1000;

%%
% Open-points (as at 10-April)
%
% # Does not converge for L_INIT > 0.59





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
L_MIN = 0.1; %mm 
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

L_INIT = 0.55;
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
% *Remark a*: Here the acceleration and deceleration phases of the muscle are
% not explicitly modeled. In that the muscle has no mass. Any acceleration
% or deceleration that arises is due to the dynamic change in the
% force-velocity relationship; i.e. $F_{VEL} = f(V,L,V_{MAX},F_{MAX}^{LOCAL})$
%
% *Remark b:* The shape parameters for the force-velocity relationship are
% computed with the $F_{MAX}^{LOCAL}$ such that;
%
% $$b = a\times V_{MAX}/F_{MAX}^{LOCAL}$$
%
% where $a=0.25$ is defined. This assures that the force-velocity
% relationship always gives $F_{MAX}^{LOCAL}$ when $V=0$.

% Specify the external load
F_EXT = 3; 
% define a time-step
tStep = 0.01; %s
tol = 0.0005;

V_C = ones(1,MEMBUFF);
L_CURR(MEMBUFF) = zeros;
F_MUSC(MEMBUFF) = zeros;
F_PASS(MEMBUFF) = zeros;
F_ACT(MEMBUFF) = zeros;
i = 1;
L_CURR(i) = L_INIT;
V_C(i) = 1;
V_PREV = 1;
L_PREV = L_MAX/2;
F_PASS(i) = F_PASS_INIT;


while abs(V_PREV) > tol && L_PREV > 0 && L_PREV < L_MAX
    % find $F_{MAX}^{LOCAL}$
    % find the "local" F_{MAX}^{LOCAL}, i.e. the F_{MAX} at the given (inital) length L_{INIT}
    F_MAXLOCAL = force_length(L_CURR(i), L_REST, F_MAX);
    F_MAXLOCAL = F_MAXLOCAL*F_MAX;

    % contractile velocity for a given F_{EXT}
    % remove the passive component
    F_EXT_ACC = F_EXT - F_PASS(i);
    V_C(i) = contractile_vel(F_EXT_ACC, V_MAX, F_MAXECC, F_MAXLOCAL);
    % update muscle length
    L_CURR(i+1) = L_CURR(i) + V_C(i)*tStep;
    V_PREV = V_C(i);
    L_PREV = L_CURR(i); 
    i = i+1;
    % find the total muscle force (without RESIDUAL FORCE ENHACENEMENT)
    [F_MUSC(i),F_ACT(i),F_PASS(i)] = force_muscle(L_CURR(i), ... 
        L_REST, V_MAX, V_C(i), F_MAXLOCAL, alpha, F_MAXECC);    
    
end

xvec = 0:tStep:tStep*(i-2);
yvec = L_CURR(1:i-1);
ftitle = 'Length of Muscle';
xtitle = 't (s)';
ytitle = 'Length (mm)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

fnum = fnum + 1;
xvec = 0:tStep:tStep*(i-2);
yvec = V_C(1:i-1);
ftitle = 'Velocity of Muscle';
xtitle = 't (s)';
ytitle = 'Velocity (mm/s)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)

fnum = fnum + 1;
xvec = L_CURR(1:i-1);
yvec = F_MUSC(1:i-1);
ftitle = 'Total Muscle Force';
xtitle = 'L (mm)';
ytitle = 'F_{MUSC} (N)';
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)
xvec = L_CURR(1:i-1);
yvec = F_ACT(1:i-1);
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)
xvec = L_CURR(1:i-1);
yvec = F_PASS(1:i-1);
plotxy(xvec, yvec, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)
legend('F_{total}','F_{active}','F_{passive}')
