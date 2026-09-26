% RUN_COMPARISON  Run the multibody and continuum CVT chain models and compare them.
%
% Purpose
%   Solves one operating point with MultibodyModel (discrete pins) and CMMModel
%   (continuum) and plots normal force, sliding angle and chain tension side by
%   side (Fig. 5 of the reference paper). The point is to see where the two
%   modelling routes agree and where they diverge.
%
% Usage
%   >> run_comparison
%   Run from this folder so that CMMModel.m and MultibodyModel.m are on the path.
%
% Inputs
%   None. The operating point is set in the params struct below.
%
% Outputs
%   One figure, 2 x 3 panels: P (N), psi (rad), F (N) against time (s);
%   odd panels = multibody, even panels = CMM. Nothing is saved to disk.
%
% Dependencies
%   CMMModel.m, MultibodyModel.m. Base MATLAB; sgtitle needs R2018b or later.
%
% Known limitations
%   - Several fields (F0, T_load, S_DR, S_DN, J, TDN, c, m, i) are set but not
%     read by either model; kp is assigned twice with the same value.
%   - sigma = 200 kg/m gives a centrifugal tension sigma*omega^2*R^2 = 1.72 kN,
%     almost equal to the 1.75 kN initial tension hard-coded in both models, so
%     F - sigma*omega^2*R^2 is only ~27 N. A metal CVT chain is of order 1 kg/m.
%     Both models' pressure levels depend strongly on this value.
%   - The two models use different time grids (50 vs 51 points), so the
%     curves are compared visually, not point by point.
%   - See the headers of CMMModel.m and MultibodyModel.m for model-level issues.

function run_comparison()
    %% Operating point
    % clear inside a function only clears this function's (empty) workspace.
    clear;clc;
    params.R      = 58.7e-3;   % chain pitch radius on the pulley (m)
    params.omega  = 50;        % pulley angular speed (rad/s); ~480 rpm
    params.sigma  = 200;       % chain mass per unit length (kg/m); see header, high for a real chain
    params.n      = 78;        % number of pins in the loop (-)
    params.m      = 0.01;      % mass of one pin (kg)
    params.k      = 5.16e4;    % longitudinal stiffness of one link (N/m)
    params.kp     = 8.63e7;    % axial stiffness of a pin against the sheaves (N/m)
    params.Ll     = 9e-3;      % pin pitch / link length (m)
    params.i      = 168e-3;    % pulley centre distance (m)
    params.c      = 0;         % link damping (N*s/m); 0 = undamped
    params.mu0    = 0.09;      % saturated friction coefficient (-); 0.07-0.12 for steel in CVT fluid
    params.vs0    = 0.001;     % friction regularisation speed (m/s); 1 mm/s keeps the law close
                               %   to Coulomb while staying smooth at vs = 0
    params.kp     = 8.63e7;    % duplicate of the kp line above (N/m)
    params.F0     = 1000;      % initial tension (N); not used, both models hard-code 1750 N
    params.T_load = 50;        % load torque (N*m); not used
    params.S_DR   = 20e3;      % driver pulley clamping force (N); not used
    params.S_DN   = 20e3;      % driven pulley clamping force (N); not used
    params.beta0  = 11*pi/180; % undeformed sheave half-angle (rad); 11 deg, typical 9-11 deg
    params.J      = 0.039;     % pulley moment of inertia (kg*m^2); not used
    params.TDN    = 50;        % driven pulley torque (N*m); not used

    %% Build and run both models on the same parameters
    cmm = CMMModel(params);
     mb = MultibodyModel(params);

    cmm_results = cmm.simulate(params);
    mb_results = mb.simulate(params);

    %% Unpack results
    % CMM: time history of one tracked chain element.
    theta = cmm_results.theta;
    normal_force_cmm = cmm_results.p;
    sliding_angle_cmm = cmm_results.psi;
    chain_force_cmm = cmm_results.F;
    te=cmm_results.t;
    % Multibody: pin states (n x N) plus the scalar histories that are plotted.
    t = mb_results.t;
    x = mb_results.x;
    y = mb_results.y;
    vx = mb_results.vx;
    vy = mb_results.vy;
    P=mb_results.P;
    F=mb_results.F;
    psi=mb_results.psi;
    % Earlier post-processing, kept for reference. It derived "force" and
    % "angle" from position and velocity magnitudes, which have the wrong
    % units; the models now return P, psi and F directly.
    % normal_force_mb = zeros(1, length(t));
    % sliding_angle_mb = zeros(1, length(t));
    % chain_force_mb = zeros(1, length(t));
    %
    % for i = 1:length(t)
    %     normal_force_mb(i) = sqrt(sum(x(:,i).^2 + y(:,i).^2));
    %     sliding_angle_mb(i) = atan2(mean(vy(:,i)), mean(vx(:,i)));
    %     chain_force_mb(i) = sqrt(mean(vx(:,i).^2 + vy(:,i).^2));
    % end

    %% Plot: multibody (blue) vs CMM (red)
    figure('Position', [100 100 1200 800]);

    % Normal force P (N)
    subplot(2,3,1);
    plot(t, P, 'b-');
    title('MB Model: Normal Force');
    xlabel('t[s]'); ylabel('P[N]');

    subplot(2,3,2);
    plot(te, normal_force_cmm, 'r-');
    title('CMM Model: Normal Force');
    xlabel('t'); ylabel('P[N]');

    % Sliding angle psi (rad)
    subplot(2,3,3);
    plot(t, psi, 'b-');
    title('MB Model: Sliding Angle');
    xlabel('t[s]'); ylabel('ψ[rad]');

    subplot(2,3,4);
    plot(te, sliding_angle_cmm, 'r-');
    title('CMM Model: Sliding Angle');
    xlabel('t'); ylabel('ψ[rad]');

    % Chain tension F (N)
    subplot(2,3,5);
    plot(t, F, 'b-');
    title('MB Model: Chain Force');
    xlabel('t[s]'); ylabel('F[N]');

    subplot(2,3,6);
    plot(te, chain_force_cmm, 'r-');
    title('CMM Model: Chain Force');
    xlabel('t'); ylabel('F[N]');

    sgtitle('Comparison of MB and CMM Models');
end
