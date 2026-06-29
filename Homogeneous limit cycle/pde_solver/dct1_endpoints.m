function Y = dct1_endpoints(X, mode)
%DCT1_ENDPOINTS  DCT-I-style transform for endpoint-including grids.
%
%   Y = dct1_endpoints(X, 'forward') converts nodal values on an
%   endpoint-including grid to cosine coefficients.
%
%   Y = dct1_endpoints(X, 'inverse') reconstructs nodal values from those
%   coefficients.
%
%   The transform acts along the second dimension. Thus X may be a row
%   vector, a column vector, or an m x N array whose rows are transformed
%   independently.
%
%   This is not MATLAB's default DCT-II normalization. It is the DCT-I-style
%   normalization used by diffusion_step_DCT.m for homogeneous Neumann
%   boundary conditions on an endpoint grid.
    if nargin < 2
        mode = 'forward';
    end

    if ~isnumeric(X) || ~isreal(X)
        error('dct1_endpoints: X must be a real numeric array.');
    end

    wasVector = isvector(X);
    wasColumn = iscolumn(X);

    if wasVector
        X = X(:).';   % work internally with row vectors
    end

    [~, N] = size(X);

    if N == 0
        error('dct1_endpoints: empty input is not allowed.');
    end

    if N == 1
        Y = X;   % trivial one-point case
        if wasVector && wasColumn
            Y = Y.';
        end
        return;
    end

    switch lower(mode)
        case {'forward', 'f'}
            Yraw = raw_dct1(X);
            Y = Yraw / (N - 1);
            Y(:, [1, N]) = 0.5 * Y(:, [1, N]);

        case {'inverse', 'i'}
            Yraw = (N - 1) * X;
            Yraw(:, [1, N]) = 2 * Yraw(:, [1, N]);
            Y = raw_dct1(Yraw) / (2 * (N - 1));

        otherwise
            error('dct1_endpoints: mode must be ''forward'' or ''inverse''.');
    end

    if wasVector && wasColumn
        Y = Y.';
    end
end

function Yraw = raw_dct1(X)
%RAW_DCT1  Unnormalized DCT-I along the 2nd dimension.
%
%   For each row x of length N,
%     Yraw_k = x_0 + (-1)^k x_{N-1} + 2*sum_{j=1}^{N-2} x_j cos(pi*j*k/(N-1)).

    [~, N] = size(X);

    % Even extension: [x0, ..., x_{N-1}, x_{N-2}, ..., x1]
    Xext = [X, X(:, N-1:-1:2)];

    F = fft(Xext, [], 2);
    Yraw = real(F(:, 1:N));
end