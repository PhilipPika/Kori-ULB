function [u_water, v_water, flux_x, flux_y, ResTime, vel_mag_eff] = fCalcWaterVel(ctr, par, flw, Wd, kegdsx, kegdsy, MASK)
% FCALCWATERVEL Calculates fluxes, velocities and residence times.
% 
% Note: To maintain mass conservation in AdvDiffTracer.m, use the
% flux_x and flux_y outputs as inputs for the tracer module.

% --- 1. Effective Depth & Gradient ---
Wd_eff = max(Wd, par.Wdmin);
grad_mag = sqrt(kegdsx.^2 + kegdsy.^2);
grad_mag_eff = max(grad_mag, 1e-9);

% --- 2. Calculate Staggered Fluxes [m2/a] ---
% [GA] These are the inputs needed for AdvDiffTracer.m
flux_x = flw .* (kegdsx ./ grad_mag_eff);
flux_y = flw .* (kegdsy ./ grad_mag_eff);

% --- 3. Calculate Velocities [m/a] ---
u_water = flux_x ./ Wd_eff;
v_water = flux_y ./ Wd_eff;

% Masking
u_water(MASK==0) = 0;
v_water(MASK==0) = 0;
flux_x(MASK==0) = 0;
flux_y(MASK==0) = 0;

% --- 4. Calculate Local Residence Time --- See Obsidian note
vel_mag = sqrt(u_water.^2 + v_water.^2);
vel_mag_eff = max(vel_mag, 1e-6); 

% Residence time: (dx * 1000) / speed [m / (m/a) = a]
ResTime = (ctr.delta * 1e3) ./ vel_mag_eff;

% Set residence time to NaN or a very high value outside the mask
ResTime(MASK == 0) = NaN;
end