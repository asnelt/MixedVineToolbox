% Copyright (C) 2016 Arno Onken
%
% This file is part of the Mixed Vine Toolbox.
%
% The Mixed Vine Toolbox is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses/>.

%% Construct 4D mixed copula vine
disp('Constructing mixed copula vine...');
d = 4; % Dimension
vine.type = 'c-vine'; % Canonical vine type
% Set margins
vine.margins = cell(d,1);
% Standard normal margin
vine.margins{1}.dist = 'norm';
vine.margins{1}.theta = [0;1];
vine.margins{1}.iscont = true; % Continuous margin
% Gamma margin
vine.margins{2}.dist = 'gam';
vine.margins{2}.theta = [2;4];
vine.margins{2}.iscont = true; % Continuous margin
% Poisson margin
vine.margins{3}.dist = 'poiss';
vine.margins{3}.theta = 10;
vine.margins{3}.iscont = false; % Discrete margin
% Binomial margin
vine.margins{4}.dist = 'bino';
vine.margins{4}.theta = [20;0.4];
vine.margins{4}.iscont = false; % Discrete margin
% Set copula families
vine.families = cell(d);
vine.theta = cell(d);
% Gaussian copula family
vine.families{1,2} = 'gaussian';
vine.theta{1,2} = 0.5;
% Student copula family
vine.families{1,3} = 'student';
vine.theta{1,3} = [0.5;2];
% Clayton copula family 
vine.families{1,4} = 'clayton';
vine.theta{1,4} = 5;
% Clayton copula family rotated 90° clockwise
vine.families{2,3} = 'claytonrot090';
vine.theta{2,3} = 5;
% Clayton copula family survival transformed
vine.families{2,4} = 'claytonrot180';
vine.theta{2,4} = 5;
% Independence
vine.families{3,4} = 'ind';
vine.theta{3,4} = [];
fprintf('\n');

%% Test probability density function
disp('Calculating probability density function on a grid...');
% Calculate probability density function on lattice
x1gv = linspace(-3,3,100);
x2gv = linspace(0.5,25,100);
x3gv = 0:20;
x4gv = 0:15;
[x1,x2,x3,x4] = ndgrid(x1gv,x2gv,x3gv,x4gv);
p = mixedvinepdf(vine,[x1(:),x2(:),x3(:),x4(:)]);
p = reshape(p,size(x1));
% Plot 2D margins
figure('Name','2D margins of mixed vine copula PDF','Position',[0,0,1600,1000]);
subplot(2,3,1);
margin12 = reshape(sum(sum(p,4),3),[length(x1gv),length(x2gv)]);
imagesc(x1gv,x2gv,margin12');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
margin13 = reshape(sum(sum(p,4),2),[length(x1gv),length(x3gv)]);
imagesc(x1gv,x3gv,margin13');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
margin14 = reshape(sum(sum(p,3),2),[length(x1gv),length(x4gv)]);
imagesc(x1gv,x4gv,margin14');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
margin23 = reshape(sum(sum(p,4),1),[length(x2gv),length(x3gv)]);
imagesc(x2gv,x3gv,margin23');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
margin24 = reshape(sum(sum(p,3),1),[length(x2gv),length(x4gv)]);
imagesc(x2gv,x4gv,margin24');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
margin34 = reshape(sum(sum(p,2),1),[length(x3gv),length(x4gv)]);
imagesc(x3gv,x4gv,margin34');
colormap('hot');
set(gca,'YDir','normal');
xlabel('Margin 3');
ylabel('Margin 4');
fprintf('\n');

%% Test sampling
disp('Sampling from mixed copula vine...');
% Draw samples
cases = 1000;
x = mixedvinernd(vine,cases);
% Plot samples in 2D
figure('Name','Mixed vine copula samples in 2D','Position',[0,0,1600,1000]);
subplot(2,3,1);
scatter(x(:,1),x(:,2),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
subplot(2,3,2);
scatter(x(:,1),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
subplot(2,3,3);
scatter(x(:,1),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 4');
subplot(2,3,4);
scatter(x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
subplot(2,3,5);
scatter(x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 4');
subplot(2,3,6);
scatter(x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 3');
ylabel('Margin 4');
% Plot samples in 3D
figure('Name','Mixed vine copula samples in 3D','Position',[0,0,1000,1000]);
subplot(2,2,1);
scatter3(x(:,1),x(:,2),x(:,3),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 3');
subplot(2,2,2);
scatter3(x(:,1),x(:,2),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 2');
zlabel('Margin 4');
subplot(2,2,3);
scatter3(x(:,1),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 1');
ylabel('Margin 3');
zlabel('Margin 4');
subplot(2,2,4);
scatter3(x(:,2),x(:,3),x(:,4),20,[0 0 0],'filled');
xlabel('Margin 2');
ylabel('Margin 3');
zlabel('Margin 4');
fprintf('\n');

%% Test vine fit
disp('Fitting parameters to samples...');
% Construct argument to specify which margins are continuous
iscont = false(d,1);
for i = 1:d
    iscont(i) = vine.margins{i}.iscont;
end
vineest = mixedvinefit(x,vine.type,iscont);
% Compare ground-truth and estimated mixed vine
fprintf('\nGround-truth mixed vine:\n');
for i = 1:d
    fprintf(' Margin %d:     %s with parameters\t',i,vine.margins{i}.dist);
    for j = 1:length(vine.margins{i}.theta)
        fprintf('\t%.2f ',vine.margins{i}.theta(j));
    end
    fprintf('\n');
end
for i = 1:d
    for j = (i+1):d
        fprintf(' Copula (%d,%d): %s with parameters\t',i,j,vine.families{i,j});
        for k = 1:length(vine.theta{i,j})
            fprintf('\t%.2f ',vine.theta{i,j}(k));
        end
        fprintf('\n');
    end
end
fprintf('\n');
disp('Estimated mixed vine:');
for i = 1:d
    fprintf(' Margin %d:     %s with parameters\t',i,vineest.margins{i}.dist);
    for j = 1:length(vineest.margins{i}.theta)
        fprintf('\t%.2f ',vineest.margins{i}.theta(j));
    end
    fprintf('\n');
end
for i = 1:d
    for j = (i+1):d
        fprintf(' Copula (%d,%d): %s with parameters\t',i,j,vineest.families{i,j});
        for k = 1:length(vineest.theta{i,j})
            fprintf('\t%.2f ',vineest.theta{i,j}(k));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Estimate entropy
disp('Estimating entropy of mixed copula vine...');
alpha = 0.05; % Significance level of estimate (95% confidence)
erreps = 1e-1; % Maximum standard error
[h,stderr] = mixedvineentropy(vine,alpha,erreps);
% Display results
disp([' Estimated entropy:          ' num2str(h) ' bit']);
disp([' Standard error of estimate: ' num2str(stderr) ' bit']);
fprintf('\n');

