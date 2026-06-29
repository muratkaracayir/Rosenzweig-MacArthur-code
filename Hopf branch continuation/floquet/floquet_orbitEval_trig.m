function evalOrbit = floquet_orbitEval_trig(UV, T, cfg)
%FLOQUET_ORBITEVAL_TRIG  Trigonometric time evaluator for one periodic orbit.
%
%   evalOrbit = floquet_orbitEval_trig(UV, T, cfg)
%
%   builds a function handle for evaluating the stored periodic orbit at
%   arbitrary times. The orbit is assumed to be stored on the periodic grid
%
%       t_j = (j-1)*T/Nt,  j = 1,...,Nt,
%
%   with the endpoint t = T excluded.
%
%   The returned handle can be called as
%
%       [u, v] = evalOrbit(t)
%
%   or, if time derivatives are also needed,
%
%       [u, v, u_t, v_t] = evalOrbit(t).
%
%   All returned quantities are column vectors of length Nx.

    validate_inputs(UV, T, cfg);

    Nx = cfg.periodic_orbit.Nx;
    Nt = cfg.periodic_orbit.Nt;
    Ntot = Nx * Nt;

    % Orbit blocks, stored as Nx-by-Nt arrays.
    U = reshape(UV(1:Ntot), Nx, Nt);
    V = reshape(UV(Ntot+1:2*Ntot), Nx, Nt);

    % Fourier coefficients in time. The normalization gives the convention
    % X(t) = sum_k FX(:,k) exp(i*omega_k*t).
    FU = fft(U, [], 2) / Nt;
    FV = fft(V, [], 2) / Nt;

    % Signed temporal Fourier frequencies.
    k = 0:(Nt - 1);
    half = floor(Nt / 2);
    kSigned = k;
    kSigned(k > half) = k(k > half) - Nt;
    omega = (2*pi/T) * kSigned;

    evalOrbit = @orbit_eval;

    function varargout = orbit_eval(t)
        if ~isscalar(t) || ~isfinite(t)
            error('floquet_orbitEval_trig:BadTime', ...
                'The evaluation time t must be a finite scalar.');
        end

        % Periodic wrapping. This also maps t = T back to t = 0.
        t = mod(t, T);

        E = exp(1i * omega * t);

        u = real(FU * E.');
        v = real(FV * E.');

        if nargout <= 2
            varargout = {u, v};
            return;
        end

        % Time derivatives of the trigonometric interpolant.
        Ew = (1i * omega) .* E;
        u_t = real(FU * Ew.');
        v_t = real(FV * Ew.');

        varargout = {u, v, u_t, v_t};
    end
end

function validate_inputs(UV, T, cfg)
%VALIDATE_INPUTS  Basic checks for the orbit evaluator.

    if ~isstruct(cfg) || ~isfield(cfg, 'periodic_orbit')
        error('floquet_orbitEval_trig:BadCfg', ...
            'The input cfg must contain cfg.periodic_orbit.');
    end

    if ~isfield(cfg.periodic_orbit, 'Nx') || ...
            ~isfield(cfg.periodic_orbit, 'Nt')
        error('floquet_orbitEval_trig:BadGrid', ...
            'The input cfg.periodic_orbit must contain Nx and Nt.');
    end

    Nx = cfg.periodic_orbit.Nx;
    Nt = cfg.periodic_orbit.Nt;

    if ~isscalar(Nx) || ~isscalar(Nt) || Nx <= 0 || Nt <= 0 || ...
            Nx ~= round(Nx) || Nt ~= round(Nt)
        error('floquet_orbitEval_trig:BadGrid', ...
            'Nx and Nt must be positive integers.');
    end

    if ~isscalar(T) || ~isfinite(T) || T <= 0
        error('floquet_orbitEval_trig:BadT', ...
            'The period T must be a positive finite scalar.');
    end

    if numel(UV) ~= 2 * Nx * Nt
        error('floquet_orbitEval_trig:BadSize', ...
            'UV must have length 2*Nx*Nt.');
    end
end
