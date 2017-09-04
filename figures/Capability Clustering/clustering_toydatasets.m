% This file is to demo the capability of CHDD in clustering
% Data: toy datasets


%% dataset of concentric rings
load('Dat_Ring.mat')
A = [Ring2 Ring3]';
A = A * scalem(A,'variance');

% c1 = chdc_train(A,0.1,0,1e-4,1e-2);          % no kernel
[c2,v2] = chdc_train(A,0.1,0.2,1e-4,1e-2);        % Guassian kernel


%% dataset of moons
load('Dat_moons.mat');
B = moons;
B = B * scalem(B,'variance');

% c3 = chdc_train(B,0.05,0,1e-4,1e-2);
[c4,v4] = chdc_train(B,0.1,0.2,1e-4,1e-2);


%% dataset of blobs
load('Dat_blobs.mat');
C = blobs;
C = C * scalem(C,'variance');

% c5 = chdc_train(C,0.05,0,1e-4,1e-2);
[c6,v6] = chdc_train(C,0.1,0.2,1e-4,1e-2);


%% dataset of blobs
load('Dat_Multi.mat');
D = Multi';
D = D * scalem(D,'variance');

% c7 = chdc_train(D,0.05,0,1e-4,1e-2);
[c8,v8] = chdc_train(D,0.5,0.18,1e-6,1e-2);


%% plot
figure

% indc_set1 = unique(c1);
% figure
% hold on
% for k=1:length(indc_set1)
%     h(k) = scatter(A(c1(1,:)==indc_set1(k),1),A(c1(1,:)==indc_set1(k),2));
% end
% hold off

indc_set2 = unique(c2);
subplot(1,4,1)
hold on
for k=1:length(indc_set2)
    scatter(A(c2(1,:)==indc_set2(k),1), A(c2(1,:)==indc_set2(k),2));
    
    index = logical(v2.W.data.si' .* (c2==indc_set2(k)));
    scatter(A(index,1), A(index,2));
end
hold off

% indc_set3 = unique(c3);
% figure
% hold on
% for k=1:length(indc_set3)
%     h(k) = scatter(B(c3(1,:)==indc_set3(k),1),B(c3(1,:)==indc_set3(k),2));
% end
% hold off

indc_set4 = unique(c4);
subplot(1,4,2)
hold on
for k=1:length(indc_set4)
    scatter(B(c4(1,:)==indc_set4(k),1),B(c4(1,:)==indc_set4(k),2));
    
    index = logical(v4.W.data.si' .* (c4==indc_set4(k)));
    scatter(B(index,1), B(index,2));
end
hold off

% indc_set5 = unique(c5);
% figure
% hold on
% for k=1:length(indc_set5)
%     h(k) = scatter(C(c5(1,:)==indc_set5(k),1),C(c5(1,:)==indc_set5(k),2));
% end
% hold off

indc_set6 = unique(c6);
subplot(1,4,3)
hold on
for k=1:length(indc_set6)
    scatter(C(c6(1,:)==indc_set6(k),1),C(c6(1,:)==indc_set6(k),2));
    
    index = logical(v6.W.data.si' .* (c6==indc_set6(k)));
    scatter(C(index,1), C(index,2));
end
hold off

% indc_set7 = unique(c7);
% figure
% hold on
% for k=1:length(indc_set7)
%     h(k) = scatter(D(c7(1,:)==indc_set7(k),1),D(c7(1,:)==indc_set7(k),2));
% end
% hold off

indc_set8 = unique(c8);
subplot(1,4,4)
hold on
for k=1:length(indc_set8)
    scatter(D(c8(1,:)==indc_set8(k),1),D(c8(1,:)==indc_set8(k),2));
    
    index = logical(v8.W.data.si' .* (c8==indc_set8(k)));
    scatter(D(index,1), D(index,2));
end
hold off