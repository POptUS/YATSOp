function [y, fvec] = calfun_sample(x, probspecs, probtype)
%     This is a modified version of the subroutine calfun.m
%     available at
%     http://www.mcs.anl.gov/~more/dfo/
%
%     Inputs:
%       x 	n-by-1 array 
%       probspecs is a struc containing:
%           m : positive integer (length of output from mghvec)
%           nprob : positive integer defining the number of the problem
%           trunc (optional) : positive scalar magnitude at which component
%               function values should be truncated
%	    sigma (optional) : scalar defining standard deviation of noise
% 	    seed  (optional) : scalar defining the rand/randn seed used
%       probtype is a string specifying the type of problem desired:
%           'smooth' corresponds to smooth (noise-free) problems
%           'absnormal' corresponds to stochastic Gaussian absolute noise
%           'absuniform' corresponds to stochastic uniform absolute noise
%           'abswild' corresponds to deterministic absolute noise
%           'relnormal' corresponds to stochastic Gaussian relative noise
%           'reluniform' corresponds to stochastic uniform relative noise
%           'relwild' corresponds to deterministic relative noise
%           'nondiff' corresponds to piecewise-smooth problems
%	**Note: the noise is applied independently to each component before
%		the components are squared and summed, additional variance
%		control will necessarily need to account for the value m
%
%     Outputs:
%       y : scalar function value at x
%       fvec : m-by-1 component function values at x
%
%     To store the evaluation history, additional variables are passed
%     through global variables. These may be commented out if a user
%     desires. They are:
%       nfev is a non-negative integer containing the number of function
%          evaluations done so far (nfev=0 is a good default).
%          after calling calfun_sample, nfev will be incremented by one.
%       np is a counter for the test problem number. np=1 is a good
%          default if only a single problem/run will be done.
%       fvals is a matrix containing the history of function
%          values, the entry fvals(nfev+1,np) being updated here.
%



%! Todo:
% Reference the starting point script for allowable definitions. Point to
% BENDFO. Provide problem descriptions

% global m nprob probtype fvals nfev np trunc


if size(x, 1) == 1
    x = x(:); % Ensure input is a column vector
end

nin = size(x, 1); % Problem dimension

eid = 'Input:dimensionIncompatible';
if nin ~= probspecs.n
    error(eid, 'Input x is not of size n by 1.');
end

% Generate the vector
fvec = mghvec(probspecs.m,probspecs.n,x,probspecs.nprob);

% Optional truncation:
if isfield(probspecs,'trunc') && max(abs(fvec))>probspecs.trunc
    fvec = sign(fvec).*min(abs(fvec),probspecs.trunc);
    %  display('Component function value exceeds trunc')
end

% Optional for noisy problems:
if ~isfield(probspecs,'sigma') % default for noisy problems
    probspecs.sigma = 10^-3;
end
if isfield(probspecs,'seed') % if seed is specified
    rand('seed',probspecs.seed)
end

% Calculate the function value
switch probtype
    case 'absnormal'
        z = probspecs.sigma*randn(probspecs.m,1);
        fvec = fvec + z;
        y = sum(fvec.^2);
    case 'absuniform'
        z = (probspecs.sigma*sqrt(3))*(2*rand(probspecs.m,1) - ones(probspecs.m,1));
        fvec = fvec + z;
        y = sum(fvec.^2);
    case 'abswild'
        z = 0.9*sin(100*norm(x,1))*cos(100*norm(x,inf)) + 0.1*cos(norm(x,2));
        z = z*(4*z^2 - 3);
        y = sum(fvec.^2) + z;
    case 'nondiff'
        y = sum(abs(fvec));
    case 'relnormal'
        z = probspecs.sigma*randn(probspecs.m,1);
        fvec = fvec.*(1 + z);
        y = sum(fvec.^2);
    case 'reluniform'
        z = (probspecs.sigma*sqrt(3))*(2*rand(probspecs.m,1) - ones(probspecs.m,1));
        fvec = fvec.*(1 + z);
        y = sum(fvec.^2);
    
    case 'absnormal2'
        z = probspecs.sigma*randn;
        y = sum(fvec.^2) + z;
    case 'absuniform2'
        z = (probspecs.sigma*sqrt(3))*(2*rand - 1);
        y = sum(fvec.^2) + z;
    case 'reluniform2'
        z = (probspecs.sigma*sqrt(3))*(2*rand - 1);
        y = sum(fvec.^2)*(1 + z);
    case 'relnormal2'
        z = probspecs.sigma*randn;
        y = sum(fvec.^2)*(1 + z);
    
        
    case 'relwild'
        z = 0.9*sin(100*norm(x,1))*cos(100*norm(x,inf)) + 0.1*cos(norm(x,2));
        z = z*(4*z^2 - 3);
        y = (1 + probspecs.sigma*z)*sum(fvec.^2);
    case 'smooth'
        y = sum(fvec.^2);
end

% Update the function value history
%nfev = nfev + 1;
%fvals(nfev,np) = y;
