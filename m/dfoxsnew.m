function [x, prob] = dfoxsnew(m, n, nprob)
%     This is a new (Jan 2022) version that allows nprob to go over 22!
%
%     This is a Matlab version of the subroutine dfoxs.f
%     This subroutine specifies the standard starting points for the
%     functions defined by subroutine dfovec as used in:
%
%     [x, prob] = dfoxsnew(m, n, nprob)
%
%     Inputs:
%       n       [int] positive input dimension
%       m       [int] positive fvec dimension (used for checking)
%       nprob   [int] positive problem number
%                       nprob<=22 are from DFO set in [1]
%     Outputs:
%       x       [n-by-1 dbl] array containing the standard starting
%                               point for problem nprob
%       prob    [struct] containing any known problem details:
%         name  [str] string for known (e.g., CUTEST) or other name
%         fbest [dbl] lower bound on SOS version of the problem
%         xbest [n-by-1 dbl] array at which fbest is attained
%         h     [dbl] forward difference parameter recommended at x for
%                       stated values of m and n
%         hm    [int] value of m for which the above h is recommended
%         hn    [int] value of n for which the above h is recommended
%
%     [1] Benchmarking Derivative-Free Optimization Algorithms. Jorge J.
%         More' & Stefan M. Wild. SIAM J. Optimization 20(1):172-191, 2009.
%     [2] MGH
%     [3] CUTEST problem definitions (sif) https://bitbucket.org/optrove/sif/src/master/
%     [4] OPM problem definitions (m) https://github.com/gratton7/OPM/tree/main/problems/
%
%     March 2022.
%     Argonne National Laboratory

% Todo (also see mghvec): reg test, errorhandling, timing, refs,
% Currently unused: FREURONE 209, +

x = zeros(n, 1); % Note: get rid of this if not a listed case.

% Error messages
eid = 'Input:dimensionIncompatible';
err{1} = 'Input m must be equal to n.';
err{2} = 'Input m must be at least n.';
err{3} = 'Input n must be even.';
err{4} = 'Input sqrt(n) must be an integer.';
err{5} = 'Input nthroot(n,3) must be an integer.';
err{6} = 'Input sqrt(n/2) must be an integer.';
err{7} = 'Input m must be equal to 2*n-2.';
err{8} = 'Input m must be equal to 2*(n-4).';
err{9} = 'Input n must be at least 5.';
err{10} = 'Input n must be at least 2.';
err{11} = 'Input (sqrt(1+4*n)-1)/2 must be an integer.';
err{12} = 'Input m must be equal to n+1.';
err{13} = 'Input m must be equal to 2*n.';
err{14} = 'Input n/4 must be an integer.';
err{15} = 'Input n must be at least 13.';
err{16} = 'Input (n+2)/3 must be an integer.';
err{17} = 'Input m must be equal to (5*n-8)/3.';
err{18} = 'Input m must be equal to n+2.';
err{19} = 'Input sqrt(n+1) must be an integer.';

switch nprob
    case 1 % MGH 32 linear function - full rank or rank 1.
        if m < n
            error(eid, err{2});
        end
        x = ones(n, 1);
        prob.name   = 'ARGLALE'; % CUTEST
        prob.h      = 2e-7;
        prob.hm     = 4000;
        prob.hn     = 2000;
        prob.fbest  = m - n; % at multiple, see MGH

    case 2 % MGH 33 linear function - full rank or rank 1.
        if m < n
            error(eid, err{2});
        end
        x = ones(n, 1);
        prob.name   = 'ARGLBLE'; % CUTEST
        prob.h      = 1e-7;
        prob.hm     = 4000;
        prob.hn     = 2000;
        prob.fbest  = m * (n - 1) / (2 * (2 * m + 1)); % at multiple, see MGH

    case 3 % MGH 34 linear function - full rank or rank 1.
        if m < n
            error(eid, err{2});
        end
        x = ones(n, 1);
        prob.name   = 'ARGLCLE'; % CUTEST
        prob.h      = 1e-7;
        prob.hm     = 4000;
        prob.hn     = 2000;
        prob.fbest  = (m^2 + 3 * m - 6) / (2 * (2 * m - 3)); % at multiple, see MGH

    case 4 % An extended version of MGH 1:  Rosenbrock
        if m ~= n
            error(eid, err{1});
        end
        x = -ones(n, 1);
        prob.name   = 'EXTROSNB'; % CUTEST & OPM
        prob.h      = 2e-8;
        prob.hm     = 1000;
        prob.hn     = 1000;
        prob.fbest  = 0;
        prob.xbest  = zeros(n, 1);

    case 5 % Another version of MGH 1:  Rosenbrock
        % SW calls this alternative chained rosenbrock
        if m ~= 2 * n - 2
            error(eid, err{7});
        end
        x = -ones(n, 1);
        prob.name   = 'ROSENBR'; % CUTEST & OPM
        prob.h      = 3e-8;
        prob.hm     = 1998;
        prob.hn     = 1000;
        prob.fbest  = 0;
        prob.xbest  = ones(n, 1);

    case 121 % MGH 21: Extended Rosenbrock (block version)
        if m ~= n
            error(eid, err{1});
        end
        if mod(n, 2)
            error(eid, err{3});
        end
        x(1) = -1.2;
        x(2) = 1;
        x = repmat(x(1:2), n / 2, 1);
        prob.name   = 'SROSENBR'; % CUTEST
        prob.h      = 6e-9;
        prob.hm     = 2000;
        prob.hn     = 2000;
        prob.fbest  = 0;
        prob.xbest  = ones(n, 1);

    case 15 % MGH 35  Chebyquad function.
        if m < n
            error(eid, err{2});
        end
        for k = 1:n
            x(k) = k / (n + 1);
        end
        prob.name   = 'CHEBYQAD'; % CUTEST
        prob.h      = 5e-11;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 16 % MGH 27 Brown almost-linear function.
        if m ~= n
            error(eid, err{1});
        end
        x = .5 * ones(n, 1);
        prob.name   = 'BROWNALE'; % CUTEST
        prob.h      = 7e-8;
        prob.hm     = 1000;
        prob.hn     = 1000;
        prob.fbest  = 0;
        prob.xbest  = ones(n, 1); % and maybe elsewhere

    case 19 % BDQRTIC
        if n < 5
            error(eid, err{5});
        end
        if m ~= 2 * (n - 4)
            error(eid, err{8});
        end
        x = ones(n, 1);
        prob.name   = 'BDQRTIC'; % CUTEST
        prob.h      = 2e-8;
        prob.hm     = 3992;
        prob.hn     = 2000;

    case 120 % CUBE, cubic version of Rosenbrock
        if m ~= n
            error(eid, err{1});
        end
        % Alt:       x = .5*ones(n,1);
        x = [-1.2; ones(n - 1, 1)];
        prob.name   = 'CUBE'; % CUTEST, OPM
        prob.h      = 3e-9;
        prob.hm     = 2000;
        prob.hn     = 2000;
        prob.fbest  = 0;
        prob.xbest  = ones(n, 1);

    case 21 % MANCINO, see opm note on exponent
        if m ~= n
            error(eid, err{1});
        end
        if n < 2
            error(eid, err{10});
        end
        for i = 1:n
            ss = 0;
            for j = 1:n
                ss = ss + (sqrt(i / j) * ((sin(log(sqrt(i / j))))^5 + (cos(log(sqrt(i / j))))^5));
            end
            x(i) = -8.710996e-4 * ((i - 50)^3 + ss);
        end
        prob.name   = 'MANCINO'; % CUTEST, OPM
        prob.h      = 2e-3;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 217 %  FREUROTH  extended freudenstein and roth function.
        if m ~= 2 * n - 2
            error(eid, err{7});
        end
        if mod(n, 2)
            error(eid, err{3});
        end
        x(1:2:n - 1) = 0.5;
        x(2:2:n) = -2;
        prob.name   = 'FREUROTH'; % CUTEST
        prob.h      = 6e-8;
        prob.hm     = 9998;
        prob.hn     = 5000;

    case 126  % VarTrig (disagrees with Table 3 ARGTRIG) MGH # 26
        if m ~= n
            error(eid, err{1});
        end
        x = ones(n, 1) / n; % Note: elsewhere you start at 1 or 1/(5n)
        prob.name   = 'VarTrig'; % OPM
        prob.h      = 2e-7;
        prob.hm     = 1000;
        prob.hn     = 1000;
        prob.fbest  = 0;
        prob.xbest  = zeros(n, 1);

    case 201  % Artificial turning point problem (doi:10.1137/0801016)
        if m ~= n
            error(eid, err{1});
        end
        x = ones(n, 1);
        prob.name   = 'ARTIF'; % CUTEST
        prob.h      = 7e-9;
        prob.hm     = 5000;
        prob.hn     = 5000;
        prob.fbest  = 0;
        prob.xbest  = zeros(n, 1);

    case 202  % ARWHDNE
        if m ~= 2 * n - 2
            error(eid, err{7});
        end
        x = ones(n, 1);
        prob.name   = 'ARWHDNE'; % CUTEST
        prob.h      = 2e-8;
        prob.hm     = 9998;
        prob.hn     = 5000;

    case 128  % BDVALUES (same func as 228, different x)
        if m ~= n
            error(eid, err{1});
        end
        h = 1 / (n + 1);
        t = [1:n]' * h;
        x = 1000 * (t .* (t - 1));
        prob.name   = 'BDVALUES'; % CUTEST
        prob.h      = 1e-7;
        prob.hm     = 1000;
        prob.hn     = 1000;
        % f*=0 ?

    case 130 % BROYDN3D Broyden Tridiagonal Function
        if m ~= n
            error(eid, err{1});
        end
        x = -ones(n, 1);
        prob.name   = 'BROYDN3D'; % CUTEST
        prob.h      = 7e-9;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 131 % MGH #31 %does not fully agree with pycutest BROYDNBD
        if m ~= n
            error(eid, err{1});
        end
        x = -ones(n, 1);
        prob.name   = 'BROYDNBD'; % CUTEST
        prob.h      = 1e-8;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 203 % BRATU2D. Problem 3 from More' 1990 collection
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n), 1) ~= 0
            error(eid, err{4});
        end
        x = zeros(n, 1);
        prob.name   = 'BRATU2D'; % CUTEST
        prob.h      = 7e-11;
        prob.hm     = 4900;
        prob.hn     = 4900;

    case 204 % BRATU2DT. Problem 3 from More' 1990 collection (same x0 as BRATU2D)
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n), 1) ~= 0
            error(eid, err{4});
        end
        x = zeros(n, 1);
        prob.name   = 'BRATU2DT'; % CUTEST
        prob.h      = 7e-11;
        prob.hm     = 4900;
        prob.hn     = 4900;

    case 205 % BRATU3D. 3D version of problem 3 from More' 1990 collection
        if m ~= n
            error(eid, err{1});
        end
        %! Need to revisit for octave versus matlab
        if mod(nthroot(n, 3), 1) ~= 0
            error(eid, err{5});
        end
        x = zeros(n, 1);
        prob.name   = 'BRATU3D'; % CUTEST
        prob.h      = 7e-10;
        prob.hm     = 3375;
        prob.hn     = 3375;

    case 206 % Complex 2D version (3.8) from More' 1990 collection
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n / 2), 1) ~= 0
            error(eid, err{6});
        end
        x = zeros(n, 1);
        prob.name   = 'CBRATU2D'; % CUTEST
        prob.h      = 8e-11;
        prob.hm     = 2888;
        prob.hn     = 2888;

    case 207 % CHANDHEQ: problem 4 from More' 1990 collection
        if m ~= n
            error(eid, err{1});
        end
        x = ones(n, 1);
        prob.name   = 'CHANDHEQ'; % CUTEST
        prob.h      = 1e-8;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 208 % EIGENB  OPM:eigenbls
        if m ~= n
            error(eid, err{1});
        end
        if mod((sqrt(1 + 4 * n) - 1) / 2, 1) ~= 0
            error(eid, err{11});
        end
        p = 0.5 * (sqrt(1 + 4 * n) - 1); % This is what must be integer
        Q = speye(p);
        x = [full(Q(:)); ones(p, 1)];
        prob.name   = 'EIGENB'; % CUTEST; OPM:eigenbls
        prob.h      = 9e-9;
        prob.hm     = 2550;
        prob.hn     = 2550;
        A = diag(2 * ones(p, 1)) - diag(ones(p - 1, 1), -1) - diag(ones(p - 1, 1), 1);
        [V, D] = eig(A);
        prob.xbest  = [reshape(V', p * p, 1); diag(D)];
        prob.fbest  = 0; % (1e-28 for n=2550 in Matlab)

    case 219 % EIGENA  OPM:eigenals
        if m ~= n
            error(eid, err{1});
        end
        if mod((sqrt(1 + 4 * n) - 1) / 2, 1) ~= 0
            error(eid, err{11});
        end
        p = 0.5 * (sqrt(1 + 4 * n) - 1); % This is what must be integer
        Q = speye(p);
        x = [full(Q(:)); ones(p, 1)];
        prob.name   = 'EIGENA'; % CUTEST; OPM:eigenals
        prob.h      = 1e-8;
        prob.hm     = 3660;
        prob.hn     = 3660;
        prob.xbest  = [reshape(eye(p), p * p, 1); [1:p]'];
        prob.fbest  = 0;

    case 220 % EIGENC  OPM:eigencls
        if m ~= n
            error(eid, err{1});
        end
        if mod((sqrt(1 + 4 * n) - 1) / 2, 1) ~= 0
            error(eid, err{11});
        end
        p = 0.5 * (sqrt(1 + 4 * n) - 1); % This is what must be integer
        Q = speye(p);
        x = [full(Q(:)); ones(p, 1)];
        prob.name   = 'EIGENC'; % CUTEST; OPM:eigencls
        prob.h      = 1e-8;
        prob.hm     = 2550;
        prob.hn     = 2550;
        A = diag(p:-1:1) + diag(ones(p - 1, 1), -1) + diag(ones(p - 1, 1), 1);
        [V, D] = eig(A);
        prob.xbest  = [reshape(V', p * p, 1); diag(D)];
        prob.fbest  = 0; % (1e-26 for n=3660 in Matlab)

    case 129 % INTEGREQ
        if m ~= n
            error(eid, err{1});
        end
        h = 1 / (n + 1);
        x = h * [1:n]' .* (h * [1:n]' - 1);
        prob.name   = 'INTEGREQ'; % CUTEST
        prob.h      = 2e-9;
        prob.hm     = 2550;
        prob.hn     = 2550;

    case 228 % MOREBV. Boundary value problem from MGH 28
        % Same func as 128, x/1000
        % OPM uses padded dims, so need to morebv('objf',[0; X0; 0]);
        if m ~= n
            error(eid, err{1});
        end
        h = 1 / (n + 1);
        t = [1:n]' * h;
        x = t .* (t - 1);
        % x=ones(n,1); % Gratton & Toint note this may be better for large n
        prob.name   = 'MOREBV'; % CUTEST & OPM
        prob.h      = 5e-13;
        prob.hm     = 1000;
        prob.hn     = 1000;
        % f*=0 ?

    case 210 % MSQRTA
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n), 1) ~= 0
            error(eid, err{4});
        end
        Npar = sqrt(n);
        xx = sin([1:n].^2);
        xx(Npar * 2 + 1) = 0; % note: opm did not do this
        x = (xx - 0.8 * sin([1:n].^2))';
        prob.name   = 'MSQRTA'; % CUTEST, also OPM
        prob.h      = 7e-8;
        prob.hm     = 4900;
        prob.hn     = 4900;
        prob.xbest  = sin([1:n].^2);
        prob.xbest(2 * Npar + 1) = 0;
        prob.fbest  = 0;

    case 211 % MSQRTB
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n), 1) ~= 0
            error(eid, err{4});
        end
        xx = sin([1:n].^2);
        x = (xx - 0.8 * sin([1:n].^2))';
        prob.name   = 'MSQRTB'; % CUTEST, also OPM
        prob.h      = 9e-8;
        prob.hm     = 1024;
        prob.hn     = 1024;
        prob.xbest  = sin([1:n].^2);
        prob.fbest  = 0;

    case 112 % OSCIGRNE
        % The roots of the gradient of Yurii Nesterov's "oscillating path" problem
        if m ~= n
            error(eid, err{1});
        end
        x = [-2; ones(n - 1, 1)];
        prob.name   = 'OSCIGRNE'; % CUTEST
        prob.h      = 4e-9;
        prob.hm     = 1000;
        prob.hn     = 1000;
        prob.xbest  = ones(n, 1);
        prob.fbest  = 0;

    case 123 % PENLT1NE (MGH 23) Penalty Function I, m=n+1;
        if m ~= n + 1
            error(eid, err{12});
        end
        x = [1:n]';
        prob.name   = 'PENLT1NE'; % CUTEST
        prob.h      = 1e-5;
        prob.hm     = 1001;
        prob.hn     = 1000;

    case 124 % Penalty2, m=2*n;
        % This a modification of PENLT2NE (MGH 24) Penalty Function II
        % It replaces the original 10 in the exponent by n.
        if m ~= 2 * n
            error(eid, err{13});
        end
        x = 0.5 * ones(n, 1);
        prob.name   = 'Penalty2'; % CUTEST
        prob.h      = 5e-9;
        prob.hm     = 2000;
        prob.hn     = 1000;

    case 113 % POWELLSG (m=n divisible by 4)
        % Extended Powell singular function. MGH 13
        if m ~= n
            error(eid, err{1});
        end
        if mod(n / 4, 1) ~= 0
            error(eid, err{14});
        end
        % Note! OPM uses [ -3; -1; 0; 1 ]
        x(1) = 3;
        x(2) = -1;
        x(3) = 0;
        x(4) = 1;
        x = repmat(x(1:4), n / 4, 1);
        prob.name   = 'POWELLSG'; % CUTEST
        prob.h      = 4e-8;
        prob.hm     = 2000;
        prob.hn     = 2000;

    case 218 % POWELLSE (m=n divisible by 4)
        % equations form of Extended Powell singular function.
        if m ~= n
            error(eid, err{1});
        end
        if mod(n / 4, 1) ~= 0
            error(eid, err{14});
        end
        x(1) = 3;
        x(2) = -1;
        x(3) = 0;
        x(4) = 1;
        x = repmat(x(1:4), n / 4, 1);
        prob.name   = 'POWELLSE'; % CUTEST
        prob.h      = 5e-8;
        prob.hm     = 1000;
        prob.hn     = 1000;

    case 213 % SPMSQRT  OPM:spmsqrt %n=3*npar-2; m=5*npar-6; npar>4
        if n < 13
            error(eid, err{15});
        end % npar>4
        if mod((n + 2) / 3, 1) ~= 0
            error(eid, err{16});
        end % npar = (n+2)/3
        if m ~= (5 * n - 8) / 3
            error(eid, err{17});
        end % m = 5*npar-6
        x = 0.2 * sin([1:n].^2)';
        prob.name   = 'SPMSQRT'; % CUTEST, opm
        prob.h      = 4e-8;
        prob.hm     = 1664;
        prob.hn     = 1000;

    case 125 % VARDIMNE (MGH 25)
        if m ~= n + 2
            error(eid, err{18});
        end
        x = 1 - [1:n]' / n;
        prob.name   = 'VARDIMNE'; % CUTEST
        prob.h      = 2e-8;
        prob.hm     = 1002;
        prob.hn     = 1000;
        prob.xbest  = ones(n, 1);
        prob.fbest  = 0;

    case 214 % YATP1SQ (Think this is YATP1LS)
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n + 1), 1) ~= 0
            error(eid, err{19});
        end
        Npar = sqrt(n + 1) - 1;
        x = [zeros(2 * Npar, 1); 6 * ones(Npar^2, 1)];
        prob.name   = 'YATP1SQ';
        prob.h      = 2e-7;
        prob.hm     = 2600;
        prob.hn     = 2600;

    case 215 % YATP2SQ
        if m ~= n
            error(eid, err{1});
        end
        if mod(sqrt(n + 1), 1) ~= 0
            error(eid, err{19});
        end
        Npar = sqrt(n + 1) - 1;
        x = [zeros(2 * Npar, 1); 10 * ones(Npar^2, 1)];
        prob.name   = 'YATP2SQ'; % CUTEST
        prob.h      = 1e-7;
        prob.hm     = 2600;
        prob.hn     = 2600;

    case 216 % ConnBand , m=n;  Problem 57 in
        %   *   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
        %   *   "Performance of a multifrontal scheme for partially separable
        %   *   optimization",
        if m ~= n
            error(eid, err{1});
        end
        x(1:2:end) = 1;
        x(2:2:end) = -1;
        prob.name   = 'ConnBand';
        prob.h      = 3e-9;
        prob.hm     = 1000;
        prob.hn     = 1000;
        prob.xbest  = zeros(n, 1);
        prob.fbest  = 0;

        % ! Remaining functions untested

    case 209 %  FREURONE   freudenstein and roth function.
        if m ~= 2 * n - 2
            error(eid, err{7});
        end
        if mod(n, 2)
            error(eid, err{3});
        end
        x(1:2:n - 1) = 0.5;
        x(2:2:n) = -2;
        prob.name   = 'FREURONE'; % CUTEST

    case 212 % SEMICN2U
        % problem 10 in "A collection of nonlinear model problems"
        x = zeros(n, 1); % Check this

    otherwise
        disp('Unknown problem.');
        clear x; % This causes error
end
