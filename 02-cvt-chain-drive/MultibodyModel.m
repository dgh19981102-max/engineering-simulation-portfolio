% MULTIBODYMODEL  Discrete (multibody) model of a metal chain on one CVT pulley.
%
% Purpose
%   Resolves the chain as n rocker pins. Each pin is pressed into the sheave
%   wedge through a linear axial spring kp (pin compression, eq. (18)); friction
%   follows a regularised Coulomb law (eq. (17)); neighbouring pins are joined
%   by longitudinal link springs k. The resulting normal force P, sliding angle
%   psi and chain tension F are compared with the continuum CMMModel in
%   run_comparison.m. Equation numbers refer to the reference paper.
%
% Usage
%   mb = MultibodyModel(params);  res = mb.simulate(params);
%   Normally called from run_comparison.m.
%
% Inputs (params struct; see the properties block for units)
%   R, omega, sigma, n, m, k, kp, Ll, i, c, mu0, vs0, beta0.
%   Delta is fixed at 1e-3 rad in the constructor.
%
% Outputs (struct returned by simulate)
%   t (s); x, y (m) and vx, vy (m/s) as n x N arrays; theta (rad, n x N);
%   F chain tension (N), P pin normal force (N), psi sliding angle (rad).
%
% Dependencies
%   Base MATLAB only.
%
% Known limitations
%   - Pin kinematics are prescribed (rigid rotation), not integrated from
%     Newton's equations; m and c are stored but not used.
%   - Pins advance at a hard-coded 20 rad/s in simulate, not obj.omega.
%   - The zone test reads pin 1 (theta(1,i)) for every j, and P(i), psi(i), F(i)
%     are overwritten for each pin, so only the last pin evaluated survives.
%   - calculateVelocities: missing parentheses make dx = x(i) - x(i-1)/theta(i)
%     - theta(i-1), and i == 1 divides by theta(1), which is 0 for pin 1 at t = 0.
%     Even with parentheses it is dx/dtheta, not a velocity.
%   - calculateChainForce receives pin j's time history for x but the whole
%     y matrix, so it indexes time rather than neighbouring links.
%   - The tension-to-pressure update uses a hard-coded friction 0.09 instead of
%     calculateFriction, and calculateFriction is fed the pitch speed
%     omega*R (~2.9 m/s) as sliding speed, so it always returns ~mu0.

classdef MultibodyModel < handle
    properties
        R;          % chain pitch radius (m)
        omega;      % pulley angular speed (rad/s)
        sigma;      % chain mass per unit length (kg/m)
        n;          % number of pins in the loop (-)
        m;          % mass of one pin (kg); not used, kinematics are prescribed
        k;          % longitudinal stiffness of one link (N/m)
        kp;         % axial stiffness of a pin against the sheaves (N/m)
        Ll;         % unstretched link length, i.e. pin pitch (m)
        i;          % pulley centre distance (m); not used in this single-pulley model
        c;          % link damping coefficient (N*s/m); not used
        mu0;        % saturated friction coefficient (-)
        vs0;        % friction regularisation speed (m/s)
        beta0;      % undeformed sheave half-angle (rad)
        Delta;      % amplitude of the sheave-angle deformation (rad)
    end

    methods
        %% Constructor
        function obj = MultibodyModel(params)
            obj.R     = params.R;
            obj.omega = params.omega;
            obj.sigma = params.sigma;
            obj.n     = params.n;
            obj.m     = params.m;
            obj.k     = params.k;
            obj.kp    = params.kp;
            obj.Ll    = params.Ll;
            obj.i     = params.i;
            obj.c     = params.c;
            obj.mu0   = params.mu0;
            obj.vs0   = params.vs0;
            obj.beta0 = params.beta0;
            obj.Delta = 0.001;  % rad; same value as CMMModel so both models see the same sheave deformation
        end

        %% Friction law, eq. (17)
        % mu = mu0*(1 - exp(-vs/vs0))   (-)
        % Regularised Coulomb friction: mu rises smoothly from 0 at vs = 0 to mu0,
        % which avoids the discontinuity of sign(vs) that makes the equations
        % stiff. vs0 sets the width of the transition; small vs0 approaches
        % pure Coulomb.
        function mu = calculateFriction(obj, vs)
            mu = obj.mu0 * (1 - exp(-vs/obj.vs0));
        end

        %% Pin contact forces, eqs. (22)-(23)
        % Tangential and radial components (N) of the two sheave contact forces
        % on a pin loaded with normal force P (N):
        %   Fh = -2*mu*P*sin(psi) / sqrt(1 + tan(beta)^2*cos(psi)^2)
        %   Fr =  2*P*sin(beta) - 2*mu*P*cos(psi) / sqrt(1 + tan(beta)^2*cos(psi)^2)
        % The factor 2 counts both sheaves. Friction opposes sliding, whose
        % direction in the pulley plane is psi (rad); the square root projects
        % the friction force from the inclined sheave face into that plane.
        function [Fh, Fr] = calculatePinForces(obj, P, psi, beta)
            beta_s = atan(tan(beta) * cos(psi));  % sheave angle along the sliding direction, eq. (9) (rad)
            denom  = sqrt(1 + tan(beta)^2 * cos(psi)^2);
            mu     = obj.calculateFriction(obj.omega * obj.R);

            Fh = -2 * mu * P * sin(psi) / denom;
            Fr = 2 * P * sin(beta) - 2 * mu * P * cos(psi) / denom;
        end

        %% Pin velocity in the rotating frame, eqs. (24)-(25)
        % (v)r =  x_dot*cos(theta) + y_dot*sin(theta)            (m/s)
        % (v)h = -x_dot*sin(theta) + y_dot*cos(theta) - r*omega   (m/s)
        % Projecting the absolute velocity onto the radial / tangential axes and
        % subtracting the sheave surface speed r*omega leaves the sliding
        % velocity of the pin relative to the sheave.
        % Here x, y, theta are the time histories of one pin and i is the time
        % index. NOTE the difference quotient lacks parentheses (see file header).
        function [vr, vh] = calculateVelocities(obj, x, y, theta,j,i)
            if i==1
                dx = x(1) ./ theta(1);
                dy = y(1) ./ theta(1);
            else

                dx = x(i)- x(i-1) ./ theta(i)-theta(i-1);
                dy = y(i)- y(i-1) ./ theta(i)-theta(i-1);
            end
            vr = dx .* cos(theta(i)) + dy .* sin(theta(i));

            % Local running radius r_i (m); differs from R once the pin moves radially.
            r = sqrt(x(i)^2 + y(i)^2);

            vh = -dx .* sin(theta(i)) + dy .* cos(theta(i)) - r * obj.omega;
        end

        %% Pin kinematics, eq. (15)
        % q_i = (x_i, y_i)^T. Pins are placed on the pitch circle R and move in
        % rigid rotation, so (vx, vy) = omega x q (m/s). This prescribes the
        % motion instead of solving for it.
        function [x, y, vx, vy] = calculatePinKinematics(obj, theta)
            x = obj.R * cos(theta);
            y = obj.R * sin(theta);

            vx = -obj.R * obj.omega * sin(theta);
            vy = obj.R * obj.omega * cos(theta);
        end

        %% Deformed sheave half-angle, eq. (2)
        % beta - beta0 = 0.5*Delta*sin(theta - thetac + pi/2)   (rad)
        % Same deformation shape as CMMModel.calculateGrooveAngle.
        function beta = calculateDeformedGrooveAngle(obj, theta, thetac)
            beta = obj.beta0 + 0.5 * obj.Delta * sin(theta - thetac + pi/2);
        end

        %% Time march
        function results = simulate(obj, params)
            % dt = 4 ms gives 51 steps over 0.2 s, matching the CMM run's horizon.
            dt = 0.004;     % time step (s)
            t  = 0:dt:0.2;  % time vector (s)
            N  = length(t);

            x     = zeros(obj.n, N);   % pin positions (m)
            y     = zeros(obj.n, N);
            vx    = zeros(obj.n, N);   % pin velocities (m/s)
            vy    = zeros(obj.n, N);
            F     = zeros(1, N);       % chain tension (N)
            P     = zeros(1, N);       % pin normal force (N)
            psi   = zeros(1, N);       % sliding angle (rad)
            theta = zeros(obj.n, N);   % pin angular positions (rad)

            % Same zones as CMMModel: 225 deg wrap, 15 deg entry zone.
            free_strand_start = pi;
            free_strand_end   = 7*pi/4;
            entry_zone_width  = pi/12;   % rad

            % Pins spaced evenly round the full circle.
            theta_init = linspace(0, 2*pi, obj.n+1);
            theta_init = theta_init(1:end-1);
            F(1) = 1750;   % initial tension (N), same as CMMModel
            P(1)=600;      % initial normal force (N), seeds the thetac integral
            for i = 1:N
                % NOTE pins advance at 20 rad/s here, not at obj.omega.
                theta(:,i) = theta_init + 20 * t(i);
                [x(:,i), y(:,i), vx(:,i), vy(:,i)] = obj.calculatePinKinematics(theta(:,i));

                for j = 1:obj.n
                    % NOTE the zone is decided by pin 1 for every j.
                    current_theta = theta(1,i);

                    is_free_strand = (current_theta >= free_strand_start) && ...
                                     (current_theta <= free_strand_end);
                    is_entry_zone  = (current_theta > free_strand_end) && ...
                                     (current_theta <= free_strand_end + entry_zone_width);

                    if is_free_strand
                        % Free strand: no sheave contact, so no normal force, no
                        % sliding and no friction to change the tension.
                        P(i)   = 0;
                        psi(i) = 0;
                        if i > 1
                            F(i) = F(i-1);
                        end
                        vr = 0;
                        vh = obj.omega * obj.R;   % pitch speed (m/s)

                    elseif is_entry_zone
                        % Entry zone: a pin seating in the groove; peak load
                        % amplified with the same heuristic ramp as CMMModel.
                        progress          = (current_theta - free_strand_end) / entry_zone_width;   % 0 -> 1
                        entry_enhancement = 1.5;   % heuristic peak amplification (-)

                        [vr, vh] = obj.calculateVelocities(x(j,:), y(j,:), theta(j,:), j, i);

                        % thetac = direction of the resultant normal force, i.e.
                        % where the wedge is pushed open most (rad).
                        num = trapz(theta(j,:), P.*sin(theta(j,:)));
                        den = trapz(theta(j,:), P.*cos(theta(j,:)));
                        thetac = atan2(num, den);
                        beta = obj.calculateDeformedGrooveAngle(current_theta, thetac);

                        % Normal force = pin spring kp times compression; a
                        % negative compression means the pin has lifted off.
                        dbi = obj.calculatePinCompression(x(j,i), y(j,i), beta);
                        if dbi >= 0
                            base_P = obj.kp * dbi;
                            P(i) = base_P * (1 + entry_enhancement * sin(pi * progress));
                        else
                            P(i) = 0;
                        end

                        psi(i) = atan2(vh, vr) * (1 + 0.5 * entry_enhancement * sin(pi * progress));

                        % Fh, Fr are computed but not used further.
                        [Fh, Fr] = obj.calculatePinForces(P(i), psi(i), beta);
                        if i > 1
                            dF = calculateChainForce(obj, x(j,:), y, i);
                            F(i) = norm(dF);
                            % Centrifugal tension sigma*omega^2*R^2 (N) and the
                            % CMM radial equilibrium; 0.09 = friction coefficient.
                            centrifugal = obj.sigma * obj.omega^2 * obj.R^2;
                            P(i) = (F(i)  - centrifugal)/(2*obj.R*(sin(obj.beta0) - ...
                            0.09*cos(obj.beta0)*cos(psi(i))));
                        end

                    else
                        % Contact arc
                        [vr, vh] = obj.calculateVelocities(x(j,:), y(j,:), theta(j,:), j, i);
                        psi(i) = pi+atan2(vh, vr);

                        num = trapz(theta(j,:), P.*sin(theta(j,:)));
                        den = trapz(theta(j,:), P.*cos(theta(j,:)));
                        thetac = atan2(num, den);
                        beta = obj.calculateDeformedGrooveAngle(current_theta, thetac);

                        % Compression is evaluated for pin 1, not pin j.
                        dbi = obj.calculatePinCompression(x(1,i), y(1,i), beta);
                        if dbi >= 0
                            P(i) = obj.kp * dbi;
                        else
                            P(i) = 0;
                        end

                        [Fh, Fr] = obj.calculatePinForces(P(i), psi(i), beta);
                        if i > 1
                            dF = calculateChainForce(obj, x(j,:), y, i);
                            F(i) = norm(dF);
                            centrifugal = obj.sigma * obj.omega^2 * obj.R^2;
                            P(i) = (F(i)  - centrifugal)/(2*obj.R*(sin(obj.beta0) - ...
                            0.09*cos(obj.beta0)*cos(psi(i))));

                        end

                    end

                end
                % Contact can push but not pull: clamp negative normal force.
                if P(i)<0
                    P(i)=0;
                end
            end

            results.t     = t;
            results.x     = x;
            results.y     = y;
            results.vx    = vx;
            results.vy    = vy;
            results.F     = F;
            results.P     = P;
            results.psi   = psi;
            results.theta = theta;
        end

        %% Pin compression, eq. (18)
        % r_i = R*tan(beta0)/tan(beta) - w_i/(2*tan(beta)) - delta_b_i/(2*tan(beta))
        % solved for the axial compression delta_b_i (m). On a wedge of
        % half-angle beta, a radial offset dr changes the sheave gap by
        % 2*tan(beta)*dr; w_i (m) is the extra gap opened by the sheave
        % deformation beta - beta0.
        function dbi = calculatePinCompression(obj, x, y, beta)
            wi=2*obj.R* tan(beta-obj.beta0);
            r = sqrt(x^2 + y^2);   % actual running radius (m)
            dbi = 2 * tan(beta) * (obj.R * tan(obj.beta0)/tan(beta) - r)-wi;
        end

        %% Link spring forces on one pin
        % Each link is a linear spring k (N/m) with rest length Ll (m); the force
        % acts along the link's unit vector tau. i < 50 keeps i+1 inside the
        % 51-step arrays. NOTE the result is the element-wise magnitude of the
        % two force vectors, not their vector sum.
        function dF = calculateChainForce(obj, x, y, i)
            if i > 1 && i<50
                dx_prev    = x(i) - x(i-1);
                dy_prev    = y(i) - y(i-1);
                alpha_prev = atan2(dy_prev, dx_prev);
                tau_prev   = [cos(alpha_prev); sin(alpha_prev)];   % unit vector of the previous link

                dx_curr    = x(i+1) - x(i);
                dy_curr    = y(i+1) - y(i);
                alpha_curr = atan2(dy_curr, dx_curr);
                tau_curr   = [cos(alpha_curr); sin(alpha_curr)];   % unit vector of the next link

                l_prev = sqrt(dx_prev^2 + dy_prev^2);   % current link lengths (m)
                l_curr = sqrt(dx_curr^2 + dy_curr^2);

                % Hooke's law on the elongation l - Ll (N).
                F_prev = -obj.k * (l_prev - obj.Ll) * tau_prev;
                F_curr = obj.k * (l_curr - obj.Ll) * tau_curr;

                dF = sqrt(F_prev.^2 + F_curr.^2);
            else
                dF=0;
            end
        end
    end
end
