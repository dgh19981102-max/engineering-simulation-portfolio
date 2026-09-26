% CMMMODEL  Continuum (CMM) model of a metal chain running on one CVT pulley.
%
% Purpose
%   Treats the chain as a continuous string of mass per length sigma wrapped on
%   a pulley whose sheaves deform elastically. Along the wrap it integrates the
%   tension ODE
%       dF/dtheta = mu*cos(beta_s)*sin(psi) / (sin(beta0) - mu*cos(beta_s)*cos(psi))
%                   * (F - sigma*omega^2*R^2)
%   and derives the sheave pressure P and sliding angle psi. Equation numbers
%   in the comments refer to the reference CMM paper. The results are compared
%   against the discrete MultibodyModel in run_comparison.m.
%
% Usage
%   cmm = CMMModel(params);  res = cmm.simulate(params);
%   Normally called from run_comparison.m.
%
% Inputs
%   params.R      pitch radius of the chain on the pulley (m)
%   params.omega  pulley angular speed (rad/s)
%   params.sigma  chain mass per unit length (kg/m)
%   beta0, Delta and mu are fixed property defaults, not read from params.
%
% Outputs (struct returned by simulate)
%   t (s), theta (rad), vr (m/s), vh (see limitations), psi (rad),
%   F chain tension (N), p pressure per unit wrap angle (N/rad).
%
% Dependencies
%   Base MATLAB only.
%
% Known limitations
%   - Kinematic march, not a solved boundary-value problem: theta = omega*t
%     follows one chain element, wrapped to [0, 2*pi).
%   - calculateTangentialVelocity returns sum(theta.*vr) over the whole arrays
%     instead of integrating eq. (1), so vh is a scalar with no clear unit.
%   - thetac (centre of the wedge opening) is only updated in the contact branch;
%     elsewhere it stays 0. trapz runs over the full, partly unfilled arrays and
%     a wrapped (non-monotonic) theta.
%   - The entry-zone amplification (1.5, sine ramp) is a heuristic shape, not
%     part of the CMM equations.
%   - The uniform deformation term of eq. (3) is omitted, and R-dot is taken as 0
%     (steady ratio).
%   - Initial tension 1750 N is hard-coded; params.F0 is not used.

classdef CMMModel < handle
    properties
        beta0 = 11 * pi/180;  % undeformed sheave half-angle (rad); 11 deg is typical of chain CVTs (9-11 deg)
        Delta = 0.001;        % amplitude of the sheave-angle deformation (rad); order 1e-3 rad
                              %   for elastic sheave bending
        mu    = 0.09;         % pin-sheave friction coefficient (-); 0.07-0.12 for steel in CVT fluid
        R;                    % chain pitch radius (m)
        omega;                % pulley angular speed (rad/s)
        sigma;                % chain mass per unit length (kg/m)
    end

    methods
        %% Constructor
        function obj = CMMModel(params)
            obj.R     = params.R;
            obj.omega = params.omega;
            obj.sigma = params.sigma;
        end

        %% Radial velocity, eq. (4)
        % vr = R_dot + a*Delta*omega*R*sin(theta - thetac)   (m/s)
        % The sheaves open where they are pushed apart most, so the chain moves
        % radially in and out once per revolution; R_dot = 0 (steady ratio).
        % a (-) is the geometric factor mapping a change of sheave angle to a
        % change of running radius on a wedge of half-angle beta0.
        function vr = calculateRadialVelocity(obj, R, theta, thetac)
            a = (1 + cos(obj.beta0)^2)/sin(2*obj.beta0);
            vr = a*obj.Delta*obj.omega*obj.R*sin(theta - thetac);
        end

        %% Tangential velocity, eq. (1)
        % Continuity of the inextensible chain: vr + d(vh)/d(theta) = 0, so vh
        % should be -integral(vr dtheta) + const. NOTE the code returns
        % sum(theta.*vr) over the whole arrays instead, a single number that does
        % not depend on the current position.
        function vh = calculateTangentialVelocity(obj, vr, theta)
            vh = sum(theta.*vr);
        end

        %% Deformed sheave half-angle, eq. (3)
        % beta - beta0 = Gamma + 0.5*Delta*sin(theta - thetac + pi/2)   (rad)
        % Sheave bending makes the groove wider on one side of the wrap, peaking
        % at thetac. The uniform term Gamma is omitted here.
        function beta = calculateGrooveAngle(obj, theta, thetac)
            beta = obj.beta0 + 0.5*obj.Delta * sin(theta - thetac + pi/2);
        end

        %% Time march
        function [results] = simulate(obj, params)
            % 50 steps over 0.2 s; at omega = 50 rad/s that is 10 rad, about 1.6
            % revolutions, wrapped to [0, 2*pi).
            t     = linspace(0, 0.2, 50);   % time (s)
            theta = params.omega * t;       % angular position of the tracked element (rad)
            theta = mod(theta, 2*pi);

            F      = zeros(size(theta));    % chain tension (N)
            P      = zeros(size(theta));    % sheave pressure per unit wrap angle (N/rad)
            vr     = zeros(size(theta));    % radial sliding velocity (m/s)
            vh     = zeros(size(theta));    % tangential velocity
            psi    = zeros(size(theta));    % sliding angle (rad)
            beta   = zeros(size(theta));    % deformed sheave half-angle (rad)
            thetac = zeros(size(theta));    % centre of the wedge opening (rad)

            % Free strand: the part of the loop between the pulleys, where the
            % chain touches no sheave. pi..7*pi/4 gives a 225 deg wrap arc.
            free_strand_start = pi;
            free_strand_end   = 7*pi/4;

            % Entry zone: the first pi/12 (15 deg) after the free strand, where
            % links enter the groove and pick up load.
            entry_zone_width = pi/12;

            F(1)      = 1750;   % tension at the start of the wrap (N)
            thetac(1) = pi/2;   % initial guess for the wedge-opening centre (rad)

            % Centrifugal tension sigma*omega^2*R^2 (N): a rotating string carries
            % this much tension with no sheave contact, so only F minus this
            % pushes the chain into the groove.
            centrifugal = obj.sigma * obj.omega^2 * obj.R^2;

            for i = 1:length(theta)
                is_free_strand = (theta(i) >= free_strand_start) && (theta(i) <= free_strand_end);
                is_entry_zone  = (theta(i) > free_strand_end) && ...
                                 (theta(i) <= free_strand_end + entry_zone_width);

                if i > 1
                    if is_free_strand
                        % Free strand: no contact, so no friction to change the tension
                        F(i)   = F(i-1);
                        P(i)   = 0;
                        psi(i) = 0;
                        vr(i)  = 0;
                        vh(i)  = obj.omega * obj.R;   % rigid transport at pitch speed (m/s)

                    elseif is_entry_zone
                        % Entry zone: contact ramps up as links seat in the groove
                        % progress runs 0 -> 1 across the zone.
                        progress = (theta(i) - free_strand_end) / entry_zone_width;

                        % Heuristic peak amplification (-): with sin(pi*progress)
                        % the pressure peaks at 1 + 1.5 = 2.5x mid-zone, mimicking
                        % the seating impact that the continuum model smears out.
                        entry_enhancement = 1.5;

                        [base_P, base_psi] = calculateBasicValues(obj, theta(i), thetac(i), F(i-1), centrifugal);

                        P(i)   = base_P * (1 + entry_enhancement * sin(pi * progress));
                        psi(i) = base_psi * (1 + 0.5 * entry_enhancement * sin(pi * progress));

                        vr(i) = obj.calculateRadialVelocity(obj.R, theta(i), thetac(i));
                        vh(i) = obj.calculateTangentialVelocity(vr, theta);
                        F(i)  = F(i-1) + calculateTensionChange(obj, theta(i), F(i-1), psi(i), centrifugal);

                    else
                        % Contact arc: full CMM equations
                        % thetac = direction of the resultant sheave pressure,
                        % i.e. where the wedge is pushed open most.
                        num = trapz(theta, P.*sin(theta));
                        den = trapz(theta, P.*cos(theta));
                        thetac(i) = atan2(num, den);

                        vr(i)   = obj.calculateRadialVelocity(obj.R, theta(i), thetac(i));
                        vh(i)   = obj.calculateTangentialVelocity(vr, theta);
                        % Sliding angle (rad): direction of the pin's sliding
                        % velocity relative to the sheave; friction acts opposite.
                        psi(i)  = pi+atan2(vh(i), vr(i));
                        beta(i) = obj.calculateGrooveAngle(theta(i), thetac(i));

                        % Sheave angle seen along the sliding direction, eq. (9) (rad).
                        beta_s = atan(tan(beta(i)) * cos(psi(i)));
                        dtheta = theta(i) - theta(i-1);
                        % Tension ODE, explicit Euler step. Only the tangential
                        % friction component (sin(psi)) changes the tension; the
                        % denominator is the radial balance of wedge normal force
                        % and the radial friction component.
                        RHS = (obj.mu * cos(beta_s) * sin(psi(i)))/(sin(obj.beta0) - ...
                            obj.mu * cos(beta_s) * cos(psi(i)));
                        F(i) = F(i-1) + dtheta * RHS * (F(i-1) - centrifugal);

                        % Radial equilibrium of an element on a wedge of half-angle
                        % beta0: 2*R*P*(sin(beta0) - mu*cos(beta0)*cos(psi)) = F - sigma*omega^2*R^2.
                        P(i) = (F(i) - centrifugal)/(2*obj.R*(sin(obj.beta0) - ...
                            obj.mu*cos(obj.beta0)*cos(psi(i))));
                    end
                end
            end

            results.theta = theta;
            results.vr    = vr;
            results.vh    = vh;
            results.psi   = psi;
            results.F     = F;
            results.p     = P;
            results.t     = t;
        end

        %% Helper: pressure and sliding angle without an updated thetac
        % Used in the entry zone. vh is fixed at pitch speed omega*R (m/s), so
        % psi is set by the radial motion alone. Same radial equilibrium as the
        % contact arc.
        function [P, psi] = calculateBasicValues(obj, theta, thetac, F, centrifugal)
            vr  = obj.calculateRadialVelocity(obj.R, theta, thetac);
            vh  = obj.omega * obj.R;
            psi = pi+atan2(vh, vr);
            P   = (F - centrifugal)/(2*obj.R*(sin(obj.beta0) - ...
                obj.mu*cos(obj.beta0)*cos(psi)));
        end

        %% Helper: tension increment in the entry zone (N)
        % Same ODE as the contact arc, but with thetac = 0 and a fixed step
        % pi/50 rad instead of the actual theta increment.
        function dF = calculateTensionChange(obj, theta, F, psi, centrifugal)
            beta   = obj.calculateGrooveAngle(theta, 0);
            beta_s = atan(tan(beta) * cos(psi));
            dtheta = pi/50;   % angular step (rad)
            RHS    = (obj.mu * cos(beta_s) * sin(psi))/(sin(obj.beta0) - ...
                obj.mu * cos(beta_s) * cos(psi));
            dF = dtheta * RHS * (F - centrifugal);
        end

    end
end
