%% Inverse involute function
% Solves tan(ang) - ang = invval for the pressure angle ang (rad) by Newton's
% method. There is no closed-form inverse. d/dx(tan x - x) = tan(x)^2, which
% vanishes at 0, so the start is kept away from 0 and each step is clamped to
% [1e-4, 1.5] rad (1.5 rad stays below the pole of tan at pi/2).
% x0: starting guess (rad); the reference transverse pressure angle is close.
function ang = solve_inv(invval, x0)
    ang = max(x0, 0.01);
    for k = 1:60
        f  = tan(ang) - ang - invval;
        df = tan(ang)^2;
        if abs(df) < 1e-12, break; end
        step = f / df;
        ang  = ang - step;
        ang  = min(max(ang, 1e-4), 1.5);
        if abs(step) < 1e-12, break; end
    end
end
