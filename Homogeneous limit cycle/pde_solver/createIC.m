function [U0, icInfo] = createIC(x, cfg, run, icIndex)
%CREATEIC  Initial conditions for homogeneous limit-cycle PDE runs.
%
%   [U0, icInfo] = createIC(x, cfg, run, icIndex)
%
%   This helper builds positive spatial initial data for the PDE time
%   stepper. The amplitudes are scaled from the supplied ODE reference
%   orbit in run.refOrbit. The default suite contains three initial
%   conditions; setting run.exploreIC = true expands this to eleven samples.
%
%   Required input:
%     x             endpoint grid, stored as a vector
%     cfg           case-parameter structure; must contain ell or L
%     run.refOrbit  reference-orbit structure with fields u and v
%     icIndex       index in the selected IC suite
%
%   Default suite, run.exploreIC = false:
%     1  equilibriumPerturbation
%     2  smoothRandom
%     3  largeSmooth
%
%   Expanded suite, run.exploreIC = true:
%     1       equilibriumPerturbation
%     2--6    smoothRandom samples
%     7--11   largeSmooth samples
    if nargin < 4
        icIndex = 1;
    end

    if ~isfield(run, 'exploreIC'), run.exploreIC = false; end
    if ~isfield(run, 'seed'),      run.seed      = 1;     end
    if ~isfield(run, 'refOrbit')
        error('createIC: run.refOrbit is required for this workflow.');
    end

    if ~isfield(cfg, 'L')
        cfg.L = cfg.ell * pi;
    end

    x = x(:).';
    L = cfg.L;

    % ------------------------------------------------------------
    % Reference-orbit scales
    % ------------------------------------------------------------
    ref = run.refOrbit;

    if ~all(isfield(ref, {'u','v'}))
        error('createIC: run.refOrbit must contain fields u and v.');
    end

    uRef = ref.u(:).';
    vRef = ref.v(:).';

    % Remove duplicated endpoint if present
    if numel(uRef) >= 2 && numel(vRef) >= 2
        if abs(uRef(end) - uRef(1)) <= 1e-12 * max(1, abs(uRef(1))) && ...
           abs(vRef(end) - vRef(1)) <= 1e-12 * max(1, abs(vRef(1)))
            uRef = uRef(1:end-1);
            vRef = vRef(1:end-1);
        end
    end

    if isempty(uRef) || isempty(vRef)
        error('createIC: reference-orbit arrays must be nonempty.');
    end

    uRefMin = min(uRef);
    uRefMax = max(uRef);
    vRefMin = min(vRef);
    vRefMax = max(vRef);

    uRefMid = 0.5 * (uRefMin + uRefMax);
    vRefMid = 0.5 * (vRefMin + vRefMax);

    uFloor = max(1e-10, 0.25 * uRefMin);
    vFloor = max(1e-10, 0.25 * vRefMin);

    uCapSmooth = max(1.5 * uRefMax, uFloor);
    vCapSmooth = max(1.5 * vRefMax, vFloor);

    uCapLarge  = max(3.0 * uRefMax, uFloor);
    vCapLarge  = max(3.0 * vRefMax, vFloor);

    % ------------------------------------------------------------
    % Select IC family and sample id
    % ------------------------------------------------------------
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
            error('With run.exploreIC=true, icIndex must be between 1 and 11.');
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
                error('With run.exploreIC=false, icIndex must be 1, 2, or 3.');
        end
    end

    % Reproducible randomness
    rng(run.seed + icIndex - 1, 'twister');

    % ------------------------------------------------------------
    % Build IC
    % ------------------------------------------------------------
    switch mode
        case 'equilibriumPerturbation'
            aU = 0.05;
            aV = 0.05;

            u0 = uRefMid * (1 + aU * cos(pi * x / L));
            v0 = vRefMid * (1 - aV * cos(pi * x / L));

            u0 = min(max(u0, uFloor), uCapSmooth);
            v0 = min(max(v0, vFloor), vCapSmooth);

        case 'smoothRandom'
            n1 = randi([1, 3]);
            n2 = randi([4, 6]);

            a1u = randSign() * (0.15 + 0.20 * rand());
            a2u = randSign() * (0.05 + 0.10 * rand());
            a1v = randSign() * (0.15 + 0.20 * rand());
            a2v = randSign() * (0.05 + 0.10 * rand());

            logPertU = a1u * cos(n1 * pi * x / L) + a2u * cos(n2 * pi * x / L);
            logPertV = a1v * cos(n1 * pi * x / L) + a2v * cos(n2 * pi * x / L);

            u0 = uRefMid * exp(logPertU);
            v0 = vRefMid * exp(logPertV);

            u0 = min(max(u0, uFloor), uCapSmooth);
            v0 = min(max(v0, vFloor), vCapSmooth);

        case 'largeSmooth'
            n1 = randi([1, 3]);
            n2 = randi([4, 6]);

            a1u = randSign() * (0.60 + 0.50 * rand());
            a2u = randSign() * (0.20 + 0.35 * rand());
            a1v = randSign() * (0.60 + 0.50 * rand());
            a2v = randSign() * (0.20 + 0.35 * rand());

            logPertU = a1u * cos(n1 * pi * x / L) + a2u * cos(n2 * pi * x / L);
            logPertV = a1v * cos(n1 * pi * x / L) + a2v * cos(n2 * pi * x / L);

            u0 = uRefMid * exp(logPertU);
            v0 = vRefMid * exp(logPertV);

            u0 = min(max(u0, uFloor), uCapLarge);
            v0 = min(max(v0, vFloor), vCapLarge);

        otherwise
            error('Unknown mode: %s', mode);
    end

    U0 = [u0; v0];

    % ------------------------------------------------------------
    % Metadata
    % ------------------------------------------------------------
    icInfo = struct();
    icInfo.mode   = mode;
    icInfo.sample = sample;

    switch mode
        case 'equilibriumPerturbation'
            icInfo.label = 'equilibriumPerturbation';

        case 'smoothRandom'
            icInfo.label = sprintf('smoothRandom-%d', sample);

        case 'largeSmooth'
            icInfo.label = sprintf('largeSmooth-%d', sample);
    end
end

% ============================================================
function s = randSign()
%RANDSIGN  Return either +1 or -1 with equal probability.
    if rand() < 0.5
        s = -1;
    else
        s = 1;
    end
end