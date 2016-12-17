%%  Single Photon Imaging Algorithm
%   Daniel Ayoub
%   ECE 436
%   Final Poject
%%  
clc; clear; close all;
fileName = 'uos-imaging/data_mannequin_face';
load(fileName);

%%
% Histogram bin settings
h_start = 2000;
h_end = 6000;
h_length = 5;
hbins = h_start:h_length:h_end;
m = length(hbins);

% Get image params
[rows, cols] = size(arrivalTimes);

% Bin data
% bindata = [];
% for i=1:rows
%     for j=1:cols
%         data = arrivalTimes{i,j};
%         data(data<h_start) = [];
%         data(data>h_end) = [];
%         bindata = [bindata; data];
%     end
% end


rms_pulsewidth = 45;
rms_pulselength_cm = rms_pulsewidth*(8e-12)*(3e8)*100;
sigs = rms_pulsewidth/h_length;
f = @(x) exp(-abs(x).^2/(2*sigs^2)); % peak is 1
% Generate dictionary
S = zeros(m,m);
t = 1:1:m;
for i=1:m        
    s = f(t-t(i));
    s = s/max(s);
    S(:,i) = s';
end

A = [S, ones(m, 1)];
delta = 1e-4;
max_iter = 10;

%% Initialize parallel pool & variables
APar = parallel.pool.Constant(A);

%%
depth = zeros(rows, cols);
tic
for i=1:10
    parfor j=1:cols
        data = arrivalTimes{i, j};
        [y, inds] = hist(data(1:15), hbins);
        sol = opt_uos(y, APar.Value, delta);
        depth(i, j) = hbins(find(sol(1:m)));
    end
end
val = toc
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