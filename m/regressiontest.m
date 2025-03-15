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

nrows = size(Var, 1);
probtype = 'smooth';
probspecs.trunc = 10^16; % Chosen so that starting point unaffected

fileID = fopen('dfof.txt', 'w+');
fprintf('num     prob     n     m            f0     hopt    time\n');
for i = 1:40 % nrows
    probspecs.nprob = Var(i, 1);
    probspecs.n = Var(i, 2);
    probspecs.m = Var(i, 3);
    factor = 10^(Var(i, 4)); % revisit!
    [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
    namestr{i} = prob.name;
    X0 = factor * X0;

    fvec = mghvec(probspecs.m, probspecs.n, X0, probspecs.nprob);
    fprintf(fileID, '%s\n', mat2str(fvec));
end
fclose(fileID);


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

fileID = fopen('dfomidf.txt', 'w+');
fprintf('num     prob     n     m            f0     hopt    time\n');
for i = 1:40 % nrows
    probspecs.nprob = Var(i, 1);
    probspecs.n = Var(i, 2);
    probspecs.m = Var(i, 3);
    factor = 10^(Var(i, 4)); % revisit!
    [X0, prob] = dfoxsnew(probspecs.m, probspecs.n, probspecs.nprob); % starting point
    namestr{i} = prob.name;
    X0 = factor * X0;

    fvec = mghvec(probspecs.m, probspecs.n, X0, probspecs.nprob);
    fprintf(fileID, '%s\n', mat2str(fvec));
end
fclose(fileID);

