%MNIST8 70000 normalised digits given by 8x8 pixels (features) in 10 classes
%
%	A = MNIST8
%
% Load the dataset in A. These are the original MNIST digits, normalised
% such that they fit exactly in images of 8x8 pixels. The first 60000 are
% the original training set, the last 10000 are the original test set.
%
% REFERENCE 
% <a href="http://yann.lecun.com/exdb/mnist/" The MNIST website</a>
%
% See also DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function a = mnist8

prdatasets(mfilename,1);
a = pr_dataset('mnist8');
a = setname(a,'MNIST8');

