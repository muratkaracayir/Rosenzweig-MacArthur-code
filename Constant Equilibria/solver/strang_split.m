function Unew = strang_split(U, cfg, Nx, dt)
%STRANG_SPLIT  One Strang-splitting step for the semidiscrete RM system.
%
%   Unew = strang_split(U, cfg, Nx, dt)
%
%   Advances the state U by one full time step dt using
%
%       diffusion(dt/2) -> reaction(dt) -> diffusion(dt/2).
%
%   The diffusion substeps are exact cosine-spectral Neumann steps, and the
%   reaction substep uses the two-stage positivity-oriented MPRK update.
%
%   Dependencies in this folder
%     diffusion_step_DCT, reaction_step_MPRK2.

    half_dt = dt / 2;

    Uhalf  = diffusion_step_DCT(U, cfg, Nx, half_dt);
    Ureact = reaction_step_MPRK2(Uhalf, cfg, dt);
    Unew   = diffusion_step_DCT(Ureact, cfg, Nx, half_dt);
end
