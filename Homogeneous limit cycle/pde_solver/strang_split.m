function Unew = strang_split(U, cfg, Nx, dt)
%STRANG_SPLIT  One Strang-splitting step for the semidiscrete PDE system.
%
%   Unew = strang_split(U, cfg, Nx, dt)
%
%   The state U is a 2 x Nx array of endpoint-grid nodal values. One full
%   step is computed as
%
%       diffusion half-step  ->  reaction step  ->  diffusion half-step.
%
%   The diffusion substeps are exact in the cosine basis for homogeneous
%   Neumann boundary conditions. The reaction substep uses the MPRK2 update
%   implemented in reaction_step_MPRK2.m.
    half_dt = dt / 2;

    Uhalf  = diffusion_step_DCT(U, cfg, Nx, half_dt);
    Ureact = reaction_step_MPRK2(Uhalf, cfg, dt);
    Unew   = diffusion_step_DCT(Ureact, cfg, Nx, half_dt);
end