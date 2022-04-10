% Sample calling syntax for the dfo and mgh funs from

clear all; % careful!

% Define Var array that shows problem specifications (nprob, n, m, xstart)
% Largescale (Table 3 from arxiv:2102.12016)
Var = [1 2000 4000 0 % ARGLALE
    2    2000 4000 0 % ARGLBLE
    201  5000 5000 0 % ARTIF
    202  5000 9998 0 % ARWHDNE
    128  1000 1000 0 % BDVALUES + Off in 5th digit but calling OK
    203  4900 4900 0 % BRATU2D
    204  4900 4900 0 % BRATU2DT
    205  3375 3375 0 % BRATU3D
    16   1000 1000 0 % BROWNALE
    130  1000 1000 0 % BROYDN3D
    131  5000 5000 0 % BROYDNBD + OK disagrees with table but agrees with OPM
    206  2888 2888 0 % CBRATU2D
    207  1000 1000 0 % CHANDHEQ
    208  2550 2550 0 % EIGENB
    217  5000 9998 0 % FREUROTH + (replaces FREURONE)
    129  1000 1000 0 % INTEGREQ
    228  1000 1000 0 % MOREBV   + (replaces MOREBVNE)
    210  4900 4900 0 % MSQRTA   + Off in 5th digit but calling OK
    211  1024 1024 0 % MSQRTB   ! Off in 3rd digit
    112  1000 1000 0 % OSCIGRNE - matches, but should confirm another value with LR
    123  1000 1001 0 % PENLT1NE
    218  1000 1000 0 % POWELLSE
    213  1000 1664 0 % SPMSQRT
    125  1000 1002 0 % VARDIMNE
    214  2600 2600 0 % YATP1SQ
    215  2600 2600 0 % YATP2SQ
    126  1000 1000 0 % VarTrig  (replaces ARGTRIG)
    216  1000 1000 0 % ConnBand (new)
    113  2000 2000 0 % POWELLSG (new)
    15   1000 1000 0 % CHEBYQAD (new)
    124  1000 2000 0 % Penalty2 (new)
    3    2000 4000 0 % ARGLCLE  (new)
    4    1000 1000 0 % EXTROSNB (new)
    5    1000 1998 0 % ROSENBR  (new)
    121  2000 2000 0 % SROSENBR (new)
    19   2000 3992 0 % BDQRTIC  (new)
    120  2000 2000 0 % CUBE     (new)
    21   1000 1000 0 % MANCINO  (new)
    219  3660 3660 0 % EIGENA   (new)
    220  2550 2550 0 % EIGENC   (new)
    212  1000 1000 0 % SEMICN2U * More' - missing
    126  1000 1000 0 % ARGTRIG  X Agrees with OPM, replaced by VarTrig
    209  5000 9998 0 % FREURONE X Disagrees, replaced by FREUROTH
    228  1000 1000 0 % MOREBVNE X Not yet in agreement, replaced with MOREBV
    ];

noiseflag = 0;

nrows = size(Var, 1);
probtype = 'smooth';
probspecs.trunc = 10^16; % Chosen so that starting point unaffected
if noiseflag
    fprintf('num     prob     n     m            f0     noise     hopt    time\n');
else
    fprintf('num     prob     n     m            f0     hopt    time\n');
end
for i = 1:40 % nrows
    probspecs.nprob = Var(i, 1);
    probspecs.n = Var(i, 2);
    probspecs.m = Var(i, 3);
    factor = 10^(Var(i, 4)); % revisit!
    [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
    namestr{i} = prob.name;
    X0 = factor * X0;

    tic;
    %    for kk=1:1e2
    y = calfun_sample(X0, probspecs, probtype);
    %   end
    ti = toc;
    % for now, overwrite y and let fvals do the work

    if noiseflag
        h = 1e-8;
        nf = 13;
        rand('state', 1); % Matlab may warn, but here's how I get reproducibility

        p = rand(probspecs.n, 1);  % The direction along which to compute derivative
        for j = 1:nf
            fval(j) = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype);
        end
        % Compute noise estimate
        [fnoise(i, 1), level, inform(i, 1)] = ECnoise(nf, fval);
        [fder2(i), s2n(i)] = f2est(@calfun_sample, nf, X0, h, p, fval, fnoise(i), probspecs, probtype);
        hopt(i, 1) = 1.68 * sqrt(fnoise(i) / abs(fder2(i)));

        p = rand(probspecs.n, 1);  % The direction along which to compute derivative
        for j = 1:nf
            fval(j) = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype);
        end
        % Compute noise estimate
        [fnoise2(i, 1), level, inform(i, 1)] = ECnoise(nf, fval);
        [fder2(i), s2n(i)] = f2est(@calfun_sample, nf, X0, h, p, fval, fnoise2(i), probspecs, probtype);
        hopt2(i, 1) = 1.68 * sqrt(fnoise2(i) / abs(fder2(i)));

        p = rand(probspecs.n, 1);  % The direction along which to compute derivative
        for j = 1:nf
            fval(j) = calfun_sample(X0 + h * (j - 7) * p, probspecs, probtype);
        end
        % Compute noise estimate
        [fnoise3(i, 1), level, inform(i, 1)] = ECnoise(nf, fval);
        [fder2(i), s2n(i)] = f2est(@calfun_sample, nf, X0, h, p, fval, fnoise3(i), probspecs, probtype);
        hopt3(i, 1) = 1.68 * sqrt(fnoise3(i) / abs(fder2(i)));

        fnoise0(i, 1) = (fnoise(i) + fnoise2(i) + fnoise3(i)) / 3;
        hopt0(i, 1) = (hopt(i) + hopt2(i) + hopt3(i)) / 3;

        fprintf('%3i  %8s  %i  %i  %12.7g %8.0e %8.0e %7.6f\n', probspecs.nprob, ...
            namestr{i}, probspecs.n, probspecs.m, y, fnoise0(i), hopt0(i), ti);
    else
        fprintf('%3i  %8s  %i  %i  %12.7g %8.0e %7.6f\n', probspecs.nprob, ...
            namestr{i}, probspecs.n, probspecs.m, y, prob.h, ti);
    end
end

if noiseflag
    % If either of these are less than 0.1, then the three estimates differ by
    % more than an order of magnitude:
    min(min([hopt hopt2 hopt3], [], 2) ./ max([hopt hopt2 hopt3], [], 2));
    min(min([fnoise fnoise2 fnoise3], [], 2) ./ max([fnoise fnoise2 fnoise3], [], 2));

    hopt = round(hopt, 1, 'significant');
    save hvals hopt namestr;
end
