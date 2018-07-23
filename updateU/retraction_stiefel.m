function Y = retraction_stiefel(X, U, t,k)
if nargin<4
    k=1;
end
if nargin < 3
    t = 1.0;
end
Y = X + t*U;
for i = 1 : k
    [Q, R] = qr(Y(:, :, i), 0);
    % The instruction with R assures we are not flipping signs
    % of some columns, which should never happen in modern Matlab
    % versions but may be an issue with older versions.
    Y(:, :, i) = Q * diag(sign(sign(diag(R))+.5));
end
end