%Breast Cancer Wisconsin (Diagnostic) Data Set
%
% http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+%28Diagnostic%29
%
% A = WDBC
% Attribute Information:
% 
% 1) ID number (ignored)
% 2) Diagnosis (M = malignant, B = benign (set as target class)) 
% 3-32) 
% 
% Ten real-valued features are computed for each cell nucleus: 
% 
% a) radius (mean of distances from center to points on the perimeter) 
% b) texture (standard deviation of gray-scale values) 
% c) perimeter 
% d) area 
% e) smoothness (local variation in radius lengths) 
% f) compactness (perimeter^2 / area - 1.0) 
% g) concavity (severity of concave portions of the contour) 
% h) concave points (number of concave portions of the contour) 
% i) symmetry 
% j) fractal dimension ("coastline approximation" - 1)

function x = breastdiag

load('wdbc.mat');

x = prdataset(wdbc,wdbclabel);
x = setname(x,'Breast Cancer Wisconsin (Diagnostic) Data Set');