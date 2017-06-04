function area = auc_pr(probs, true_labels)

% Computes the area under the precision-recall curve.
% - probs is the probability of having the correct label.
% - true_labels is the true label (0 or 1).
%
% Copyright (C) 2012 Rodolphe Jenatton, Nicolas Le Roux, Antoine Bordes, Guillaume Obozinski
%         
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% prob_1 is the probability that the label is 1.
prob_1 = probs.*true_labels + (1 - probs).*(1 - true_labels);

[dummy, isort] = sort( prob_1, 'descend');
precision = cumsum( true_labels(isort) ) ./ (1:numel(probs));
recall    = cumsum( true_labels(isort) ) / sum( true_labels );
area = trapz([0,recall],[1,precision]);
