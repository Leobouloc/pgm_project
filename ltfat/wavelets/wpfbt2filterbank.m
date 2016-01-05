function [g,a] = wpfbt2filterbank( wtdef, varargin)
%WPFBT2FILTERBANK  WPFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wpfbt2filterbank(wtdef)
%
%   Input parameters:
%         wtdef : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   WPFBT2FILTERBANK(wtdef) calculates the impulse responses g and the
%   subsampling factors a of non-iterated filterbank, which is equivalent
%   to the wavelet packet filterbank tree described by wtdef. The returned
%   parameters can be used directly in FILTERBANK, UFILTERBANK or
%   FILTERBANK.
%
%   The function internally calls WFBTINIT and passes wtdef and all
%   additional parameters to it.
%
%   Examples:
%   ---------
%
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:
%
%     [g,a] = wpfbt2filterbank({'db10',3,'dwt'});
%     filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:
%
%     [g,a] = wpfbt2filterbank({'db10',3,'full'});
%     filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   See also: wfbtinit
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wpfbt2filterbank.php

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


complainif_notenoughargs(nargin,1,'WPFBT2FILTERBANK');

definput.import = {'wfbtcommon'};
definput.flags.interscaling={'intsqrt','intnoscale','intscale'};
[flags]=ltfatarghelper({},definput,varargin);

% build the tree
wt = wfbtinit({'strict',wtdef},flags.forder);

wt = comp_wpfbtscale(wt,flags.interscaling);

nIdx = nodesLevelsBForder(wt);
% Now we need to walk the tree by levels
g = {};
a = [];
for ii=1:numel(nIdx)
    rangeLoc = cellfun(@(eEl) 1:numel(eEl.h),wt.nodes(nIdx{ii}),...
                       'UniformOutput',0);
    rangeOut = cellfun(@(eEl) numel(eEl.h),wt.nodes(nIdx{ii}));
    rangeOut = mat2cell(1:sum(rangeOut),1,rangeOut);
    [gtmp,atmp] = nodesMultid(nIdx{ii},rangeLoc,rangeOut,wt);
    g(end+1:end+numel(gtmp)) = gtmp;
    a(end+1:end+numel(atmp)) = atmp;
end
g = g(:);
a = a(:);

















