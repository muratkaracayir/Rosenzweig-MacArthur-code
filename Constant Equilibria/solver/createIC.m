function [U0, icInfo] = createIC(x, cfg, run, icIndex)
%CREATEIC  Generate initial data for constant-equilibrium simulations.
%
%   [U0, icInfo] = createIC(x, cfg, run, icIndex)
%
%   This helper constructs one of the initial conditions used by the public
%   constant-equilibrium time-stepping workflow for the diffusive
%   Rosenzweig--MacArthur system.
%
%   Input
%     x        1 x Nx or Nx x 1 grid vector on [0,L]
%     cfg      case structure with fields k, theta, lambda, u_star, v_star
%              and either L or ell
%     run      run-options structure. Relevant fields are:
%                run.exploreIC  false: use 3 basic ICs; true: use 11 ICs
%                run.seed       seed for reproducible random ICs
%                run.eqType     'coexistence', 'preyonly', or auto if empty
%     icIndex  index of the requested initial condition
%
%   Initial-condition families
%     equilibriumPerturbation  small cosine perturbation of the target state
%     smoothRandom             low-mode random smooth perturbation
%     largeSmooth              larger smooth positive perturbation
%
%   If run.exploreIC is false:
%       icIndex = 1,2,3 correspond to the three families above.
%
%   If run.exploreIC is true:
%       icIndex = 1      equilibriumPerturbation
%       icIndex = 2..6   smoothRandom, five realizations
%       icIndex = 7..11  largeSmooth, five realizations
%
%   Output
%     U0       2 x Nx initial state, with U0(1,:) = u0 and U0(2,:) = v0
%     icInfo   structure describing the selected family and target state

    if nargin < 4
        icIndex = 1;
    end

    if ~isfield(run, 'exploreIC'), run.exploreIC = false; end
    if ~isfield(run, 'seed'),      run.seed      = 1;     end

    if ~isfield(run, 'eqType') || isempty(run.eqType)
        if cfg.lambda < cfg.k
            run.eqType = 'coexistence';
        else
            run.eqType = 'preyonly';
        end
    end

    if ~isfield(cfg, 'L')
        if ~isfield(cfg, 'ell')
            error('createIC: cfg must contain either cfg.L or cfg.ell.');
        end
        cfg.L = cfg.ell * pi;
    end

    x = x(:).';
    L = cfg.L;

    % --- Target equilibrium / reference predator scale ---
    switch lower(run.eqType)
        case 'coexistence'
            if cfg.lambda >= cfg.k
                error(['createIC: coexistence initial data requested, but ', ...
                       'lambda >= k. Use run.eqType = ''preyonly''.']);
            end
            ueq = cfg.u_star;
            veq = cfg.v_star;
            vscale = max(veq, 0.1 * cfg.k / max(cfg.theta, eps));

        case 'preyonly'
            ueq = cfg.k;
            veq = 0;
            vscale = max(1, cfg.k / max(cfg.theta, eps));

        otherwise
            error('createIC: unknown run.eqType: %s.', run.eqType);
    end

    % --- Select IC family and sample id ---
    if run.exploreIC
        if icIndex == 1
            mode = 'equilibriumPerturbation';
            sample = 1;
        elseif icIndex >= 2 && icIndex <= 6
            mode = 'smoothRandom';
            sample = icIndex - 1;
        elseif icIndex >= 7 && icIndex <= 11
            mode = 'largeSmooth';
            sample = icIndex - 6;
        else
            error('createIC: with run.exploreIC=true, icIndex must be 1,...,11.');
        end
    else
        switch icIndex
            case 1
                mode = 'equilibriumPerturbation';
                sample = 1;
            case 2
                mode = 'smoothRandom';
                sample = 1;
            case 3
                mode = 'largeSmooth';
                sample = 1;
            otherwise
                error('createIC: with run.exploreIC=false, icIndex must be 1, 2, or 3.');
        end
    end

    % --- Reproducible randomness ---
    rng(run.seed + icIndex - 1, 'twister');

    % --- Build IC ---
    switch mode
        case 'equilibriumPerturbation'
            u0 = ueq * (1 + 0.05 * cos(pi * x / L));

            if strcmpi(run.eqType, 'coexistence')
                v0 = veq * (1 - 0.05 * cos(pi * x / L));
            else
                v0 = 0.05 * vscale * (1 + cos(pi * x / L));
            end

        case 'smoothRandom'
            n1 = randi([1, 3]);
            n2 = randi([4, 6]);

            a1u = randSign() * (0.05 + 0.10 * rand());
            a2u = randSign() * (0.02 + 0.08 * rand());
            a1v = randSign() * (0.05 + 0.10 * rand());
            a2v = randSign() * (0.02 + 0.08 * rand());

            u0 = ueq * (1 + a1u * cos(n1 * pi * x / L) + ...
                            a2u * cos(n2 * pi * x / L));

            if strcmpi(run.eqType, 'coexistence')
                v0 = veq * (1 + a1v * cos(n1 * pi * x / L) + ...
                                a2v * cos(n2 * pi * x / L));
            else
                base = 0.10 * vscale;
                v0 = base * (1 + a1v * cos(n1 * pi * x / L) + ...
                                 a2v * cos(n2 * pi * x / L));
            end

        case 'largeSmooth'
            n1 = randi([1, 3]);
            n2 = randi([4, 6]);

            A1u = 0.30 + 0.40 * rand();
            A2u = 0.15 + 0.25 * rand();
            A1v = 0.30 + 0.40 * rand();
            A2v = 0.15 + 0.25 * rand();

            u0 = ueq * (0.25 + ...
                        A1u * (1 + cos(n1 * pi * x / L)) + ...
                        A2u * (1 + cos(n2 * pi * x / L)));

            if strcmpi(run.eqType, 'coexistence')
                v0 = veq * (0.25 + ...
                            A1v * (1 + cos(n1 * pi * x / L)) + ...
                            A2v * (1 + cos(n2 * pi * x / L)));
            else
                v0 = vscale * (0.10 + ...
                               A1v * (1 + cos(n1 * pi * x / L)) + ...
                               A2v * (1 + cos(n2 * pi * x / L)));
            end

        otherwise
            error('createIC: unknown IC mode: %s.', mode);
    end

    % --- Positivity safeguard ---
    u0 = max(u0, eps);
    v0 = max(v0, eps);

    U0 = [u0; v0];

    icInfo.mode     = mode;
    icInfo.sample   = sample;
    icInfo.eqType   = run.eqType;
    icInfo.targetEq = [ueq; veq];
    icInfo.label    = sprintf('%s_%02d', mode, sample);
end

function s = randSign()
%RANDSIGN  Return -1 or +1 with equal probability.

    if rand() < 0.5
        s = -1;
    else
        s = 1;
    end
end
