function Lc = treeOutLen(L,doNoExt,wt)
%TREESUB
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfbtmanip/treeOutLen.php

% Copyright (C) 2005-2014 Peter L. Soendergaard <soender@users.sourceforge.net>.
% This file is part of LTFAT version 1.4.5
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

% All nodes with at least one final output.
termN = find(noOfNodeOutputs(1:numel(wt.nodes),wt)~=0);
% Range in filter outputs
outRange = rangeInLocalOutputs(termN,wt);
cRange = cell2mat(cellfun(@(rEl) rEl(:),rangeInOutputs(termN,wt),...
                  'UniformOutput',0));


Lctmp = nodeOutLen(termN,L,outRange,doNoExt,wt);
Lc = zeros(size(Lctmp));

Lc(cRange) = Lctmp;








