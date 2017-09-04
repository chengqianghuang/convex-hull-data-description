% This file is to demo the capability of CHDD in one-class classification
% Data: toy datasets

%% dataset of concentric rings
load('Dat_Ring.mat')
A = [Ring2 Ring3]';
A = A * scalem(A,'variance');

W1 = chdd(A,0.1,0,1e-4);            % no kernel
W2 = chdd(A,0.1,0.3,1e-4);          % Guassian kernel


%% dataset of moons
load('Dat_moons.mat');
B = moons;
B = B * scalem(B,'variance');

W3 = chdd(B,0.1,0,1e-4);
W4 = chdd(B,0.1,0.3,1e-4);


%% dataset of blobs
load('Dat_blobs.mat');
C = blobs;
C = C * scalem(C,'variance');

W5 = chdd(C,0.1,0,1e-4);
W6 = chdd(C,0.1,0.3,1e-4);


%% dataset of blobs
load('Dat_Multi.mat');
D = Multi';
D = D * scalem(D,'variance');

W7 = chdd(D,0.1,0,1e-4);
W8 = chdd(D,0.1,0.3,1e-4);


%% plot results
figure
subplot(2,4,1)                % circles linear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W1.data.sv(:,1),W1.data.sv(:,2),'ro')
scatter(A(:,1),A(:,2),'b+')
plotc(W1)
hold off

subplot(2,4,5)                % circle nonlinear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W2.data.sv(:,1),W2.data.sv(:,2),'ro')
scatter(A(:,1),A(:,2),'b+')
plotc(W2)
hold off

subplot(2,4,2)                % moons nonlinear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W3.data.sv(:,1),W3.data.sv(:,2),'ro')
scatter(B(:,1),B(:,2),'b+')
plotc(W3)
hold off

subplot(2,4,6)                % moons nonlinear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W4.data.sv(:,1),W4.data.sv(:,2),'ro')
scatter(B(:,1),B(:,2),'b+')
plotc(W4)
hold off

subplot(2,4,3)                % blobs linear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W5.data.sv(:,1),W5.data.sv(:,2),'ro')
scatter(C(:,1),C(:,2),'b+')
plotc(W5)
hold off

subplot(2,4,7)                % blobs nonlinear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W6.data.sv(:,1),W6.data.sv(:,2),'ro')
scatter(C(:,1),C(:,2),'b+')
plotc(W6)
hold off

subplot(2,4,4)                % square linear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W7.data.sv(:,1),W7.data.sv(:,2),'ro')
scatter(D(:,1),D(:,2),'b+')
plotc(W7)
hold off

subplot(2,4,8)                % square nonlinear
hold on
axis([-2.5 2.5 -2.5 2.5])
scatter(W8.data.sv(:,1),W8.data.sv(:,2),'ro')
scatter(D(:,1),D(:,2),'b+')
plotc(W8)
hold off