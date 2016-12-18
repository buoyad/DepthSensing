%%  Single Photon Imaging Algorithm
%   Daniel Ayoub
%   ECE 436
%   Final Poject
%%  
clc; clear; close all;
fileName = 'uos-imaging/data_mannequin_face';
load(fileName);

%%
% Camera parameters
h_start = 2000;
h_end = 6000;
h_length = 5;
hbins = h_start:h_length:h_end;
m = length(hbins);
concurrent = false; % Use parallel pool to process data

% Get image params
[rows, cols] = size(arrivalTimes);


%
% Signal generation - taken
% from sample implementation
% for output signal params
%
rms_pulsewidth = 45;
sigs = rms_pulsewidth/h_length;
f = @(x) exp(-abs(x).^2/(2*sigs^2)); 
% Generate dictionary
S = zeros(m,m);
t = 1:1:m;
for i=1:m        
    s = f(t-t(i));
    s = s/max(s);
    S(:,i) = s';
end

A = [S, ones(m, 1)];
delta = 1e-3;

% Initialize parallel pool & variables
if concurrent == true
    APar = parallel.pool.Constant(A);
    hbinsPar = parallel.pool.Constant(hbins);
end
%%
depth = zeros(rows, cols);
rtstart = tic;
functime = 0;
if concurrent == false
    for i=1:1
        for j=1:cols
            data = arrivalTimes{i, j};
            [y, ~] = hist(data(1:15), hbins);
            fstart = tic;
            [sol, t] = opt_uos(y, A, delta);
            functime = functime + t;
            depth(i, j) = hbins(find(sol(1:m)));
        end
    end
else
    for i=1:10
        parfor j=1:cols
            data = arrivalTimes{i, j};
            [y, ~] = hist(data(1:15), hbins);
            sol = find(opt_uos(y, APar.Value, delta));
            depth(i, j) = hbins(sol(1));
        end
    end
end
runningTime = toc(rtstart)

%%
load([fileName '_truth']);
depth_true = cell2mat(D_true);
depth_error = abs(depth_true - depth);
subplot(1, 3, 2);
imagesc(depth,[3550, 3700]);
axis image; colorbar; colormap(lines);
title('Generated Map');
subplot(1, 3, 1);
imagesc(depth_true, [3550, 3700]);
axis image; colorbar; colormap(spring);
title('True Map');
subplot(1, 3, 3);
imagesc(depth_error);
axis image; colorbar; colormap(summer);
title('Error');