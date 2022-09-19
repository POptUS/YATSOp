function fvec = mghvec(m, n, x, nprob)
%     This subroutine supplements More' Garbow Hillstrom problems
%     as well as other cuter/cutest and variable dimension nls problems
%
%       fvec is an output array of length m that contains the nprob
%         function evaluated at x.
%       m and n are positive integer input variables. n must not
%         exceed m; allowable sizes depend on the problem.
%       x is an input with columns of length n; allowable sizes depend on the problem.
%       nprob is a positive integer input variable which defines the
%         number of the MGH problem (or other number)
%
%     Argonne National Laboratory
%     Xiaoqian Liu and Stefan Wild. January 2022.

% Initialize things:
fvec = zeros(m, 1);

switch nprob
     case 1 % MGH 32 ARGLALE Linear function - full rank.
            % m>= n
            s = 0;
        for j = 1:n
            s = s + x(j);
        end
        temp = 2 * s / m + 1;
        for i = 1:m
            fvec(i) = -temp;
            if i <= n
                fvec(i) = fvec(i) + x(i);
            end
        end
    case 2 % MGH 33 ARGLBLE Linear function - rank 1.
        % m>= n
        s = 0;
        for j = 1:n
            s = s + j * x(j);
        end
        for i = 1:m
            fvec(i) = i * s - 1;
        end
    case 3 %  MGH 34 ARGLCLE   Linear function - rank 1 with zero columns and rows.
        % m>= n
        s = 0;
        for j = 2:n - 1
            s = s + j * x(j);
        end
        for i = 1:m - 1
            fvec(i) = (i - 1) * s - 1;
        end
        fvec(m) = -1;
    case 4 % An extended version of MGH 1: EXTROSNB Rosenbrock
        % OPM: y = extrosnb('objf',x)
        % m=n
        fvec(1) = x(1);
        for i = 2:n
            fvec(i) = 10 * (x(i)^2 - x(i - 1));
        end

    case 5 % Another version of MGH 1:  Rosenbrock ROSENBR
        % OPM: y = rosenbr('objf',x)
        % m=2*n-2
        % SW calls this alternative chained rosenbrock
        for i = 1:n - 1
            fvec(i) = 10 * (x(i)^2 - x(i + 1));
            fvec(n - 1 + i) = x(i) - 1;
        end

     case 121 % MGH 21: SROSENBR Extended Rosenbrock (block version)
        % m=n, n even
        for i = 1:n / 2
            fvec(2 * i - 1) = 10 * (x(2 * i) - x(2 * i - 1)^2);
            fvec(n / 2 + i) = 1 - x(2 * i - 1);
        end

    case 15 % MGH 35 CHEBYQAD    chebyquad function.
        % m>=n
        for j = 1:n
            t1 = 1;
            t2 = 2 * x(j) - 1;
            t = 2 * t2;
            for i = 1:m
                fvec(i) = fvec(i) + t2;
                th = t * t2 - t1;
                t1 = t2;
                t2 = th;
            end
        end
        iev = -1;
        for i = 1:m
            fvec(i) = fvec(i) / n;
            if iev > 0
                fvec(i) = fvec(i) + 1 / (i^2 - 1);
            end
            iev = -iev;
        end

    case 16 % MGH 27 BROWNALE (m=n)   Brown almost-linear function.
        sum1 = -(n + 1);
        prod1 = 1;
        for j = 1:n
            sum1 = sum1 + x(j);
            prod1 = x(j) * prod1;
        end
        for i = 1:n - 1
            fvec(i) = x(i) + sum1;
        end
        fvec(n) = prod1 - 1;

   case 19 % BDQRTIC
        % n>=5, m = (n-4)*2
        for i = 1:n - 4
            fvec(i) = (-4 * x(i) + 3.0);
            fvec(n - 4 + i) = (x(i)^2 + 2 * x(i + 1)^2 + 3 * x(i + 2)^2 + ...
                4 * x(i + 3)^2 + 5 * x(n)^2);
        end

case 120 % CUBE %OPM CUBE
        % n>=2; m=n;
        fvec(1) = (x(1) - 1.0);
        for i = 2:n
                fvec(i) = 10 * (x(i) - x(i - 1)^3);
        end
    case 21 % MANCINO
        % n >=2; m=n
        for i = 1:n
            % ss=0;
            xi2 = x(i)^2;
            v2 = sqrt(xi2 + i ./ [1:n]');
            lv2 = log(v2);
            ss = sum(v2 .* ((sin(lv2)).^5 + (cos(lv2)).^5));
            % for j=1:n
            %    v2 = sqrt(xi2 +i/j);
            %    lv2 = log(v2);
            %    ss = ss+v2*((sin(lv2))^5 + (cos(lv2))^5);
            % end
            fvec(i) = 1400 * x(i) + (i - 50)^3 + ss;
        end

    case 126 % Trigonometric function, n=m. MGH # 26
        % VarTrig (disagrees with Table 3 ARGTRIG) MGH # 26
        % OPM: argtrig https://github.com/gratton7/OPM/blob/main/problems/argtrig.m
        % y=argtrig('objf',ones(1000,1))
        % Note that OPM disagrees with MGH (has negative in front of i)
        temp = n - sum(cos(x));
        for i = 1:m
            fvec(i) = temp + i * (1 - cos(x(i))) - sin(x(i));
        end
        % f*=0 at x*=0;

    case 128 % BDVALUES Discrete boundary value problem, m=n %same as 228
        % Rheinboldt "Some Nonlinear Testproblems" (prob 6) replaces ^3 by ^2
        h = 1 / (n + 1);
        for i = 2:n - 1
            fvec(i) = 2 * x(i) - x(i - 1) - x(i + 1) + 0.5 * h^2 * (x(i) + i * h + 1)^3;
        end
        fvec(1) = 2 * x(1) - x(2) + 0.5 * h^2 * (x(1) + h + 1)^3;
        fvec(n) = 2 * x(n) - x(n - 1) + 0.5 * h^2 * (x(n) + n * h + 1)^3;

    case 130 % BROYDN3D Broyden Tridiagonal Function, m=n
        % xn= [-1; x; -1]; %wrong
        xn = [0; x; 0];
        for i = 1:m
            fvec(i) = (3 - 2 * xn(i + 1)) * xn(i + 1) - xn(i) - 2 * xn(i + 2) + 1;
        end

    case 131 % Broyden Banded Function
        % MGH #31 %does not fully agree with pycutest BROYDNBD
        fvec(1) = x(1) * (2 + 5 * x(1)^2) + 1 - x(2) .* (1 + x(2));
        fvec(m) = x(m) * (2 + 5 * x(m)^2) + 1 - sum(x(n - 5:n - 1) .* (1 + x(n - 5:n - 1)));
        x2 = x .* (1 + x);
        for i = 2:5
            term = sum(x2(1:i - 1)) + x2(i + 1);
            fvec(i) = x(i) * (2 + 5 * x(i)^2) + 1 - term;
        end
        for i = 6:m - 1
            term = sum(x2(i - 5:i - 1)) + x2(i + 1);
            fvec(i) = x(i) * (2 + 5 * x(i)^2) + 1 - term;
        end
    case 201 %     Artificial turning point problem (doi:10.1137/0801016), m=n
        x = [0; x; 0];
        for i = 1:m
            fvec(i) = atan(sin(x(i + 1) * mod(i, 100))) - (x(i) + x(i + 1) + x(i + 2)) / 20;
        end
    case 202 % ARWHDNE , m=2n-2;  Problem 55 in
        %   *   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
        %   *   "Performance of a multifrontal scheme for partially separable
        %   *   optimization",
        % but with a squared term on the linear (3-4*x(i))
        % Note that OPM uses original formulation
        % OPM: https://github.com/gratton7/OPM/blob/main/problems/arwhead.m
        for i = 1:n - 1
            fvec(i) = x(i)^2 + x(n)^2;
            % fvec(n-1+i) = sqrt(3-4*x(i)); %!?OPM Warning: not for use in NLS.
            fvec(n - 1 + i) = (3 - 4 * x(i)); % modified for nonlinear equations
        end
    case 203 % BRATU2D. Problem 3 from More' 1990 collection
        % n must be a square, m = n
        d = sqrt(n);
        h = 1 / (d + 1);
        lam = 4;  % LAMBDA is the Bratu problem parameter.  It should be positive.
        lh2 = lam * h^2;
        u = zeros(d + 2, d + 2);
        u(2:d + 1, 2:d + 1) = reshape(x, d, d);
        ff = zeros(d, d);
        for i = 1:d
            for j = 1:d
                ff(i, j) = lh2 * exp(u(i + 1, j + 1)) + (u(i + 2, j + 1) + ...
                    u(i, j + 1) + u(i + 1, j + 2) + u(i + 1, j) - 4 * u(i + 1, j + 1));
            end
        end
        fvec = ff(:);
    case 204 % BRATU2DT % Same as BRATU2D but at 6.80812=lambda (the turning point)
        % n must be a square, m = n
        d = sqrt(n);
        h = 1 / (d + 1);
        lam = 6.80812;  % lambda is the Bratu problem parameter.  It should be positive.
        lh2 = lam * h^2;
        u = zeros(d + 2, d + 2);
        u(2:d + 1, 2:d + 1) = reshape(x, d, d);
        ff = zeros(d, d);
        for i = 1:d
            for j = 1:d
                ff(i, j) = lh2 * exp(u(i + 1, j + 1)) + (u(i + 2, j + 1) + ...
                    u(i, j + 1) + u(i + 1, j + 2) + u(i + 1, j) - 4 * u(i + 1, j + 1));
            end
        end
        fvec = ff(:);

    case 205 % BRATU3D. 3D version of problem 3 from More' 1990 collection; see BRATU2D
        % n must be a cuberoot, m = n
        d = round(nthroot(n, 3)); % Warning: round
        h = 1 / (d + 1);
        lam = 6.80812;  % LAMBDA is the Bratu problem parameter.  It should be positive.
        lh2 = lam * h^2;
        u = zeros(d + 2, d + 2, d + 2);
        u(2:d + 1, 2:d + 1, 2:d + 1) = reshape(x, d, d, d);
        ff = zeros(d, d, d);
        for i = 1:d
            for j = 1:d
                for k = 1:d
                    ff(i, j, k) = lh2 * exp(u(i + 1, j + 1, k + 1)) + (u(i + 2, j + 1, k + 1) + ...
                        u(i, j + 1, k + 1) + u(i + 1, j + 2, k + 1) + u(i + 1, j, k + 1) + ...
                        u(i + 1, j + 1, k + 2) + u(i + 1, j + 1, k) - 6 * u(i + 1, j + 1, k + 1));
                end
            end
        end
        fvec = ff(:);

    case 206 % CBRATU2D Complex 2D version (3.8) of problem 3 from More' 1990 collection; see BRATU2D
        % n must be half a square root, m = n
        d = sqrt(n / 2);
        h = 1 / (d + 1);
        lam = 5;  % LAMBDA is the Bratu problem parameter.  It should be positive.
        lh2 = lam * h^2;
        x = reshape(x, d, 2 * d);
        u1 = zeros(d + 2, d + 2);
        u2 = zeros(d + 2, d + 2);
        u1(2:d + 1, 2:d + 1) = x(:, 1:d);
        u2(2:d + 1, 2:d + 1) = x(:, d + 1:2 * d);
        ff1 = zeros(d, d);
        ff2 = zeros(d, d);
        for i = 1:d
            for j = 1:d
                ff1(i, j) = lh2 * exp(u1(i + 1, j + 1)) * cos(u2(i + 1, j + 1)) + ...
                    (u1(i + 2, j + 1) + u1(i, j + 1) + u1(i + 1, j + 2) + ...
                    u1(i + 1, j) - 4 * u1(i + 1, j + 1));
                ff2(i, j) = lh2 * exp(u1(i + 1, j + 1)) * sin(u2(i + 1, j + 1)) + ...
                    (u2(i + 2, j + 1) + u2(i, j + 1) + u2(i + 1, j + 2) + ...
                    u2(i + 1, j) - 4 * u2(i + 1, j + 1));
            end
        end
        fvec = [ff1(:); ff2(:)];

    case 207 % CHANDHEQ: problem 4 from More' 1990 collection; %m=n
        % n = length(x);
        c = 1;
        xx = [1:n]' / n;
        hcw = 0.5 * c / n;
        for i = 1:n
            fvec(i) = x(i) - sum((x(i) * xx(i) * hcw) * x(:) ./ (xx(i) + xx(:))) - 1;
        end

    case 208 % EIGENB  OPM:eigenbls
        % m = n = p*(p+1)
        p = round(0.5 * (sqrt(1 + 4 * n) - 1)); % This must be integer
        Q = reshape(x(1:p^2), p, p);
        D = diag(x(p^2 + 1:n));
        lowerinds = tril(true(size(D)));
        % Note: ignore upper triangular part of A
        % Case b: A tridiagonal with 2's on diagonal, -1's on off-diagonals
        A = diag(2 * ones(p, 1)) - diag(ones(p - 1, 1), -1);
        R = A - Q' * D * Q;
        fvec(1:(p + 1) * p / 2) = R(lowerinds);
        R2 = eye(p) - Q' * Q;
        fvec(1 + (p + 1) * p / 2:m) = R2(lowerinds);

    case 219 % EIGENA  OPM:eigenals
        % m = n = p*(p+1)
        p = round(0.5 * (sqrt(1 + 4 * n) - 1)); % This must be integer
        Q = reshape(x(1:p^2), p, p);
        D = diag(x(p^2 + 1:n));
        lowerinds = tril(true(size(D)));
        % Case a: A diagonal
        A = diag(1:p);
        R = A - Q' * D * Q;
        fvec(1:(p + 1) * p / 2) = R(lowerinds);
        R2 = eye(p) - Q' * Q;
        fvec(1 + (p + 1) * p / 2:m) = R2(lowerinds);

    case 220 % EIGENC  OPM:eigencls
        % m = n = p*(p+1)
        p = round(0.5 * (sqrt(1 + 4 * n) - 1)); % This must be integer
        Q = reshape(x(1:p^2), p, p);
        D = diag(x(p^2 + 1:n));
        lowerinds = tril(true(size(D)));
        % Note: ignore upper triangular part of A
        % Case c?: A tridiagonal suggested by Wilkinson
        A = diag(p:-1:1) + diag(ones(p - 1, 1), -1);
        R = A - Q' * D * Q;
        fvec(1:(p + 1) * p / 2) = R(lowerinds);
        R2 = eye(p) - Q' * Q;
        fvec(1 + (p + 1) * p / 2:m) = R2(lowerinds);

    case 228 % MOREBV. Boundary value problem from MGH 28 (same as 128)
        % OPM uses padded dims, so need to morebv('objf',[0; X0; 0]);
        h = 1 / (n + 1);
        for i = 2:n - 1
            fvec(i) = 2 * x(i) - x(i - 1) - x(i + 1) + 0.5 * h^2 * (x(i) + i * h + 1)^3;
        end
        fvec(1) = 2 * x(1) - x(2) + 0.5 * h^2 * (x(1) + h + 1)^3;
        fvec(n) = 2 * x(n) - x(n - 1) + 0.5 * h^2 * (x(n) + n * h + 1)^3;

        %        for i=2:n-1
        %            gvec(i) = fvec(i)*(2 + 1.5*h^2*(x(i)+i*h+1)^2) - fvec(i-1) - fvec(i+1);
        %        end
        %        gvec(1) = fvec(1)*(2 + 1.5*h^2*(x(1)+h+1)^2) - fvec(2);
        %        gvec(n) = fvec(n)*(2 + 1.5*h^2*(x(n)+n*h+1)^2) - fvec(n-1);
        %        fvec=gvec(:);

    case 129 % INTEGREQ. Discrete integral equation from MGH 29, m=n
        h = 1 / (n + 1);
        hi = [1:n]' * h;
        him = 1 - hi;
        y = (x + hi + 1).^3;
        t1 = zeros(n, 1);
        t2 = zeros(n, 1);
        for i = 1:m
            t1(i) = hi(1:i)' * y(1:i);
            t2(i) = him(i + 1:n)' * y(i + 1:n);
        end
        fvec = x + (h / 2) * (him .* t1 + hi .* t2);

    case 210 % MSQRTA
        % m=n, n must be a square
        Npar = sqrt(n);
        xm = reshape(x, Npar, Npar)';
        xs = reshape(sin([1:n].^2), Npar, Npar)';
        xs(3, 1) = 0;
        diff = xs * xs - xm * xm;
        fvec = diff(:);

    case 211 % MSQRTB
        Npar = sqrt(n);
        xm = reshape(x, Npar, Npar)';
        xs = reshape(sin([1:n].^2), Npar, Npar)';
        diff = xs * xs - xm * xm;
        fvec = diff(:);

    case 112 % OSCIGRNE % m=n
        rho = 500;
        % The roots of the gradient of Yurii Nesterov's "oscillating path" problem
        % https://bitbucket.org/optrove/sif/src/master/OSCIGRNE.SIF
        % Nesterov's "Chebyshev-Rosenbrock function"?
        % fvec(1) = (x(1)-1)/2;
        % for i=1:n-1
        % fvec(i+1) = sqrt(rho)*(x(i+1)-2*x(i)^2+1);
        % end

        % ! Warning: this is Jarre's On Nesterov's Smooth Chebyshev-Rosenbrock Function
        % ! with a rho/2 for the first gradient component to match Table 3!!!

        fvec(1) = (x(1) - 1) / 2 - 4 * rho * (x(2) - 2 * x(1)^2 + 1) * x(1);
        % fvec(1)=fvec(1)/2;
        for i = 2:n - 1
            fvec(i) = 2 * rho * (x(i) - 2 * x(i - 1)^2 + 1 - 4 * x(i) * (x(i + 1) - 2 * x(i)^2 + 1));
        end
        fvec(n) = 2 * rho * (x(n) - 2 * x(n - 1)^2 + 1);

    case 123 % PENLT1NE (MGH 23) Penalty Function I, m=n+1;
        ar = sqrt(1e-5);
        for i = 1:m - 1
            fvec(i) = ar * (x(i) - 1);
        end
        fvec(m) = sum(x.^2) - 0.25;

    case 124 % Penalty2, m=2*n;
        % This a modification of PENLT2NE (MGH 24) Penalty Function II
        % It replaces the original 10 in the exponent by n.
        c = n; % Originally c=10;
        ar = sqrt(1e-5);
        fvec(1) = (x(1) - 0.2);
        for i = 2:n
            temp = exp(i / c) + exp((i - 1) / c);
            fvec(i) = ar * (exp(x(i) / c) + exp(x(i - 1) / c) - temp);
            fvec(n + i - 1) = ar * (exp(x(i) / c) - exp(-1 / c));
        end
        fvec(m) = [n:-1:1] * (x.^2) / (0.01 * c^2) - 1;

    case 113 %    % POWELLSG (m=n divisible by 4)
       % Extended Powell singular function. MGH 13
        %       https://www.sfu.ca/~ssurjano/powell.html
        % Note! OPM uses [ -3; -1; 0; 1 ]
        s5 = sqrt(5);
        s10 = sqrt(10);
        for j = 0:n / 4 - 1
            j4 = 4 * j;
            fvec(j4 + 1) = x(j4 + 1) + 10 * x(j4 + 2);
            fvec(j4 + 2) = s5 * (x(j4 + 3) - x(j4 + 4));
            fvec(j4 + 3) = (x(j4 + 2) - 2 * x(j4 + 3))^2;
            fvec(j4 + 4) = s10 * (x(j4 + 1) - x(j4 + 4))^2;
        end
    case 218 % POWELLSE (m=n divisible by 4)
       % equations form of Extended Powell singular function.
        for j = 0:n / 4 - 1
            j4 = 4 * j;
            %{
            % SW: Grad
            t1 = x(j4 + 1) + 10*x(j4 + 2);
            t2 = x(j4 + 3) -    x(j4 + 4);
            t3 = x(j4 + 2) -  2*x(j4 + 3);
            t4 = x(j4 + 1) -    x(j4 + 4);
            fvec(j4 + 1) =    t1 + 20*t4^3;
            fvec(j4 + 2) = 10*t1 +  2*t3^3;
            fvec(j4 + 3) =  5*t2 -  4*t3^3;
            fvec(j4 + 4) = -5*t2 - 20*t4^3;
            %}

            fvec(j4 + 1) = x(j4 + 1) + 10 * x(j4 + 2);
            fvec(j4 + 2) = (x(j4 + 3) -    x(j4 + 4)) * 5;
            fvec(j4 + 3) = (x(j4 + 2) -  2 * x(j4 + 3))^2;
            fvec(j4 + 4) = 10 * (x(j4 + 1) -    x(j4 + 4))^2;
        end
    case 212 % SEMICN2U
        % problem 10 in "A collection of nonlinear model problems"
        % n = length(x); % (must have 0.9*n = integer)
        % n is number discretization points
        ln = 9 * n / 10; % index of the last negative discretization point

        lambda = 0.2;   % continuation parameter
        a = -0.00009;   % = t(0) interval lower bound
        b = 0.00001;    % = t(n+1) interval upper bound
        ua = 0;         % = u(a) boundary value
        ub = 700;       % = u(b) boundary value
        ca = 1e12;
        cb = 1e13;
        beta = 40;

        t = a:(b - a) / (n + 1):b; % discretization
        % can always check that v(ln+1) is negative and v(ln+2) is not.

    case 213 % SPMSQRT
        % n=3*npar-2; m=5*npar-6; npar>4
        npar = round((n + 2) / 3);

        % m=5*npar-6; % check this was input
        for cs = 1:2
            if cs == 2 % In second pass compute optimal value
                fstore = fvec;
                x = sin([1:n].^2)';
            end

            i = 1;
            j = 0;
            % compute the function value for each element.
            for item = 1:m
                kk = 13;
                % if npar==4, kk=0; end %we require npar>4
                j = j + 1;
                if item == 4 || item == 8 || item == kk
                    i = i + 1;
                    j = 1;
                end

                num10 = 1;
                while num10
                    num10 = 0;
                    if j == 1
                        jshift = i - 3;
                        kshift = jshift;
                        if jshift <= 0
                            jshift = 0;
                        end
                        j = j + jshift;
                    end
                    if j - kshift == 6 || j == npar + 1
                        i = i + 1;
                        j = 1;
                        num10 = 1;
                    end
                end
                % compute the left and right numbers of the i-th row.
                % compute the top and bottom numbers of the j-th column.
                il = 3 * (i - 1);
                % ill=il;
                ir = i * 3 - 1;
                if ir > n
                    ir = n;
                end
                jt = 3 * (j - 1) - 1;
                % jtt=jt;
                jb = 3 * j;
                if jb > n
                    jb = n;
                end
                ishift = i - 2;
                jshift = j - 2;
                if i == 1
                    il = 1;
                    % ill=il;
                    ir = 2;
                    ishift = 0;
                end
                if j == 1
                    jt = 1;
                    % jtt=jt;
                    jb = 3;
                    jshift = 0;
                end
                % compute the product of row and column vectors.
                if ishift <= jshift
                    n1 = jshift - ishift;
                    il = il + n1;
                    s = 0.0;
                    for k2 = il:ir % 20
                        s = s + x(k2) * x(jt);
                        jt = jt + 2;
                        if jt > jb
                            break % go to 30
                        end
                    end % 20
                    % 30
                    fvec(item) = s;
                else
                    n1 = ishift - jshift;
                    jt = jt + 2 * n1;
                    s = 0.0;
                    for k2 = jt:2:jb % 40
                        s = s + x(k2) * x(il);
                        il = il + 1;
                    end % 40
                    % continue
                    fvec(item) = s;
                end
                % f= f+ s1**2
            end % item do
        end
        fvec = fstore - fvec;

    case 125 % VARDIMNE (MGH 25)
        for i = 1:n
            fvec(i) = (x(i) - 1);
        end
        fvec(n + 1) = sum([1:n]' .* (x - 1));
        fvec(n + 2) = fvec(n + 1)^2;

    case 214 % YATP1SQ (Think this is YATP1LS)
        Npar = sqrt(n + 1) - 1;
        A = 10;
        y = x(1:Npar);
        z = x(Npar + 1:2 * Npar);
        xm = reshape(x(2 * Npar + 1:end), Npar, Npar);
        isum = zeros(Npar, 1);
        jsum = zeros(Npar, 1);
        fv = zeros(Npar, Npar);
        for i = 1:Npar
            ypz = y(i) + z(i);
            isum(i, 1) = sum(sin(xm(i, :)) ./ xm(i, :)) - 1;
            jsum(i, 1) = sum(sin(xm(:, i)) ./ xm(:, i)) - 1;

            for j = 1:Npar
                fv(i, j) = xm(i, j)^3 - A * xm(i, j)^2 - ...
                    ypz *  xm(i, j) * cos(xm(i, j) - sin(xm(i, j)));
            end
        end
        fvec = [fv(:); isum; jsum];

    case 215 % YATP2SQ (Think this is YATP2LS)
        Npar = sqrt(n + 1) - 1;
        A = 1;
        y = x(1:Npar);
        z = x(Npar + 1:2 * Npar);
        xm = reshape(x(2 * Npar + 1:end), Npar, Npar);
        isum = zeros(Npar, 1);
        jsum = zeros(Npar, 1);
        fv = zeros(Npar, Npar);
        for i = 1:Npar
            ypz = y(i) + z(i);
            isum(i, 1) = sum(xm(i, :) + sin(xm(i, :))) - 1;
            jsum(i, 1) = sum(xm(:, i) + sin(xm(:, i))) - 1;

            for j = 1:Npar
                fv(i, j) = xm(i, j) - ypz * (1 + cos(xm(i, j))) - A;
            end
        end
        fvec = [fv(:); isum; jsum];

    case 216 % ConnBand , m=n;  Problem 57 in
        %   *   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
        %   *   "Performance of a multifrontal scheme for partially separable
        %   *   optimization",
        for i = 1:n - 2
            fvec(i) = (x(i) + x(i + 1) + x(n))^2;
        end
        fvec(n - 1) = x(1) - x(2);
        fvec(n) = x(n - 1) - x(n);

    case 217 % FREUROTH (agrees with OPM) m = 2*n-2, n even
        % Version of extended Freudentstein and Roth
        % problem 2 in MGH
        for i = 1:n - 1
            fvec(i) = -13 + x(i) + ((5 - x(i + 1)) * x(i + 1) - 2) * x(i + 1);
            fvec(n - 1 + i)   = -29 + x(i) + ((1 + x(i + 1)) * x(i + 1) - 14) * x(i + 1);
        end
        % ! Remaining functions untested
    case 209 % (2?)%FREURONE
        % Freudentstein and Roth; nonlinear equation version of FREUROTH
        % problem 2 in MGH
        % extended version (FREUROTH?):
        % for i =1:n/2
        % fvec(2*i-1) = -13 + x(2*i-1) + ((5 - x(2*i))*x(2*i) - 2)*x(2*i);
        % fvec(2*i)   = -29 + x(2*i-1) + ((1 + x(2*i))*x(2*i) - 14)*x(2*i);
        % end

        % (agrees with OPM) m = 2*n-2
        for i = 1:n - 1
            fvec(i) = -13 + x(i) + ((5 - x(i + 1)) * x(i + 1) - 2) * x(i + 1);
            fvec(n - 1 + i)   = -29 + x(i) + ((1 + x(i + 1)) * x(i + 1) - 14) * x(i + 1);
        end

        % f1vec(1) = 2*fvec(1);
        % f1vec(n) = 2*fvec(n);
        % for i = 2:n-1
        %    f1vec(i) = 2*fvec(i) + 2*fvec(i-1)*(10*x(i)-3*x(i)^2-2);
        %    f1vec(n-1+i) = 2*fvec(n-1+i) + 2*fvec(n-2+i)*(3*x(i)^2+2*x(i)-14);
        % end
        % fvec=f1vec(:);

        %{
        f1vec = zeros(n,1);
        f2vec = zeros(n,1);
        for i =1:n/2
            num = (-13 + x(2*i-1) + ((5 - x(2*i))*x(2*i) - 2)*x(2*i));
            f1vec(2*i-1) = 2*num;
            f1vec(2*i) = 2*num*(10*x(2*i)-3*x(2*i)^2-2);

            %f1vec(2*i) = (-13 + x(2*i-1) + ((5 - x(2*i))*x(2*i) - 2)*x(2*i));
            %f2vec(i) = (-29 + x(2*i-1) + ((1 + x(2*i))*x(2*i) - 14)*x(2*i));
            num2 = (-29 + x(2*i-1) + ((1 + x(2*i))*x(2*i) - 14)*x(2*i));
            f2vec(2*i-1) = 2*num2;
            f2vec(2*i) = 2*num2*(3*x(2*i)^2+2*x(2*i)-14);
        end
        fvec=[f1vec, f2vec]'; fvec = fvec(1:end-2); % Fix for dims not, not yet working
        %}

end

if length(fvec) ~= m
    disp('Wrong m');
    n;
    m;
    length(fvec);

end
