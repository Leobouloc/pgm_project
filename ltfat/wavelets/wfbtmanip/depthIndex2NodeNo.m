function [nodeNo,nodeChildIdx] = depthIndex2NodeNo(d,k,wt)

if(d==0)
    nodeNo=0;
    nodeChildIdx=0;
    return;
end

% find ordered nodes at depth d-1
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfbtmanip/depthIndex2NodeNo.php

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
nodesNo = getNodesInDepth(d,wt);
if(isempty(nodesNo))
   error('%s: Depth of the tree is less than given d.',mfilename); 
end

% k is index in children of ordered nodes at depth d

nodeNo = zeros(numel(k),1);
nodeChildIdx = zeros(numel(k),1);
chNo = cumsum(cellfun( @(nEl) length(nEl.g),wt.nodes(nodesNo)));
chNoZ = [0;chNo(:)];

for kIdx=1:numel(k)
    ktmp = k(kIdx);
    idx = find(chNo>ktmp,1);
    if isempty(idx)
       error('%s: Index k=%i out of bounds.',mfilename,ktmp); 
    end    
    nodeNo(kIdx) = nodesNo(idx);
    nodeChildIdx(kIdx) = ktmp-chNoZ(idx)+1;
end

% ktemp = k;
% chNo = cellfun( @(nEl) numel(nEl.g),wt.nodes(nodesNo));
% for ii=1:length(nodesNo)
%     if(ktemp<chNo(ii))
%         nodeChildIdx = ktemp+1;
%         nodeNo = nodesNo(ii);
%         if mynodeNo~=nodeNo || mynodeChildIdx ~= nodeChildIdx
%             error('mas to tu spatne');
%         end
%         return;
%     else
%         ktemp = ktemp-chNo(ii);
%     end
% end
% error('%s: Index k out of bounds.',mfilename);




function nodd = getNodesInDepth(d,wt)
% find all nodes with d steps to the root ordered
if d==1
    % return root
    nodd = find(wt.parents==0);
    return;
end    

nbf = nodesBForder(wt);
nbfTmp = nbf;
tempd = 0;
while tempd<d
    nbf(nbfTmp==0) = [];
    nbfTmp(nbfTmp==0) = [];
    nbfTmp = wt.parents(nbfTmp);
    tempd = tempd+1;
end
nodd = nbf(nbfTmp==0);

% 
% nbfPar = wt.parents(nbf);
% 
% for ii=1:numel(nbf)
%    
% end

% nodd = [];
% toGoTrough = {};
% 
% 
%    nodeNo = find(wt.parents==0);
%    toGoTrough = cell(d+1,1);
%    toGoTrough{1} = nodeNo;
%    tempd = 1;
% 
% 
% while(tempd<d)
%     
%     for jj=1:length(toGoTrough{tempd})
%        actNod = toGoTrough{tempd}(jj);
%        childrenIdx = find(wt.children{actNod}~=0);
%        ch = wt.children{actNod}(childrenIdx);
%        toGoTrough{tempd+1} = [toGoTrough{tempd+1},ch];
%     end
% 
%     tempd = tempd+1;
% end
% 
% nodd=toGoTrough{d};

