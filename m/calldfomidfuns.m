% Sample calling syntax for the dfo and mgh funs from
%  https://github.com/POptUS/YATSOp
% After cloning the above repository, please set the following variable:
%yatsop_location = '../../../YATSOp'; % For all
%yatsop_location = 'C:/Users/liljo/Dropbox/Dzahini_Wild/YATSOp'; % For KJD
%yatsop_location = '~/repos/poptus/YATSOp/'; % This location is for SW
%addpath([yatsop_location, '/m/']);



% Define Var array that shows problem specifications (nprob, n, m, h)
% Largescale (Table 3 from arxiv:2102.12016)
Var = [1 100 200 1e-7 % ARGLALE
    2    100 200 5e-8 % ARGLBLE
    3    100 200 5e-8 % ARGLCLE  (new)
    201  100 100 4e-9 % ARTIF
    202  100 198 1e-8 % ARWHDNE
    19   100 192 9e-9 % BDQRTIC  (new)
    128  100 100 1e-6 % BDVALUES + Off in 5th digit but calling OK
    203  100 100 4e-10 % BRATU2D
    204  100 100 5e-10 % BRATU2DT
    205  125 125 1e-9 % BRATU3D
    16   100 100 4e-8 % BROWNALE
    130  100 100 4e-9 % BROYDN3D
    131  100 100 4e-9 % BROYDNBD + OK disagrees with table but agrees with OPM
    206  98 98 4e-10 % CBRATU2D
    207  100 100 6e-9 % CHANDHEQ
    15   100 100 4e-10 % CHEBYQAD (new)
    216  100 100 5e-9 % ConnBand (new)
    120  100 100 4e-9 % CUBE     (new)
    219  110 110 2e-8 % EIGENA   (new)
    208  110 110 6e-9 % EIGENB
    220  110 110 1e-8 % EIGENC   (new)
    4    100 100 1e-8 % EXTROSNB (new)
    217  100 198 2e-8 % FREUROTH + (replaces FREURONE)
    129  100 100 1e-9 % INTEGREQ
    21   100 100 2e-7 % MANCINO  (new)
    228  100 100 1e-11 % MOREBV   + (replaces MOREBVNE)
    210  100 100 2e-8 % MSQRTA   + Off in 5th digit but calling OK
    211  100 100 2e-8 % MSQRTB   ! Off in 3rd digit
    112  100 100 9e-9 % OSCIGRNE - matches, but should confirm another value with LR
    124  100 200 5e-9 % Penalty2 (new)
    123  100 101 8e-7 % PENLT1NE
    218  100 100 3e-8 % POWELLSE
    113  100 100 2e-8 % POWELLSG (new)
    5    100 198 1e-8 % ROSENBR  (new)
    213  100 164 3e-8 % SPMSQRT
    121  100 100 4e-9 % SROSENBR (new)
    125  100 102 1e-8 % VARDIMNE
    126  100 100 3e-7 % VarTrig  (replaces ARGTRIG)
    214  99 99 9e-8 % YATP1SQ
    215  99 99 9e-8 % YATP2SQ
    212  100 100 0 % SEMICN2U * More' - missing
    126  100 100 0 % ARGTRIG  X Agrees with OPM, replaced by VarTrig
    209  100 198 0 % FREURONE X Disagrees, replaced by FREUROTH
    228  100 100 0 % MOREBVNE X Not yet in agreement, replaced with MOREBV
    ];

noiseflag = 0;

nrows = size(Var, 1);
probtype = 'smooth';
probspecs.trunc = 10^16; % Chosen so that starting point unaffected
if noiseflag
    addpath('~/repos/randprojections21/src/testfuncs/')
    fprintf('num     prob     n     m            f0     noise     hopt    time\n');
else
    fprintf('Problem     n    m            f0     h   \n');
%        fprintf('num     prob     n     m            f0     hopt    time\n');

end

for i = 1:40 % nrows
    probspecs.nprob = Var(i, 1);
    probspecs.n = Var(i, 2);
    probspecs.m = Var(i, 3);
    factor = 10^(0); % revisit!
    [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
    namestr{i} = prob.name;
    X0 = factor * X0;
    
    tic;
    y = calfun_sample(X0, probspecs, probtype);
    ti = toc;
    % for now, overwrite y and let fvals do the work
    
    if noiseflag
        h = 1e-11;
        nf = 23;
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

        %fprintf('%3i  %8s  %i  %i  %12.7g %8.0e %8.0e %7.6f\n', probspecs.nprob, ...
        fprintf('%8s  %i  %i  %12.7g %8.0e %8.0e \n', ...
            namestr{i}, probspecs.n, probspecs.m, y, fnoise0(i), hopt0(i));
    else
        %fprintf('%3i  %8s  %i  %i  %12.7g %8.0e %7.6f\n', probspecs.nprob, ...
        fprintf('%8s  %i  %i  %12.7g %8.0e \n', ...
            namestr{i}, probspecs.n, probspecs.m, y, prob.h);
     %   load hvals_mid
       % load hvals
        %load Bestmid2
        %fprintf('%8s &  %i & %i & %12.7g & %12.7g & %8.0e   \\\\ \n', ...
        %    namestr{i}, probspecs.n, probspecs.m, y, Best{i}.fv, Var(i,4)/100);
            %    rat(i) = (Var(i,4)/100)/(max(hopt(i)/1,1e-16));
    end
end

if noiseflag
    % If either of these are less than 0.1, then the three estimates differ by
    % more than an order of magnitude:
    min(min([hopt hopt2 hopt3], [], 2) ./ max([hopt hopt2 hopt3], [], 2));
    min(min([fnoise fnoise2 fnoise3], [], 2) ./ max([fnoise fnoise2 fnoise3], [], 2));

    hopt = round(hopt, 1, 'significant');
    save hvals_mid hopt namestr;
end
