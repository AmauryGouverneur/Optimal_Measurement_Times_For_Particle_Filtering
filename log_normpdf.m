function y = log_normdpdf(x,mu,sigma)
%LOG NORMPDF Normal probability density function (pdf).
%   Y = LOG_NORMPDF(X,MU,SIGMA) returns the log of the pdf of the normal distribution with
%   mean MU and standard deviation SIGMA, evaluated at the values in X.
%   The size of Y is the common size of the input arguments.  A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.


if nargin<1
    error(message('stats:normpdf:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    y = -0.5 * ((x - mu)./sigma).^2 ;
catch
    error(message('stats:normpdf:InputSizeMismatch'));
end
