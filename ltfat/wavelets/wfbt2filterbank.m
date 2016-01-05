function [g,a] = wfbt2filterbank( wt, varargin)
%WFBT2FILTERBANK  WFBT equivalent non-iterated filterbank
%   Usage: [g,a] = wfbt2filterbank(wt)
%
%   Input parameters:
%         wt : Wavelet filter tree definition
%
%   Output parameters:
%         g   : Cell array containing filters
%         a   : Vector of sub-/upsampling factors
%
%   [g,a]=WFBT2FILTERBANK(wt) calculates the impulse responses g and the
%   subsampling factors a of non-iterated filterbank, which is equivalent
%   to the wavelet filterbank tree described by wt. The returned
%   parameters can be used directly in FILTERBANK ant other routines.
%
%   [g,a]=WFBT2FILTERBANK({w,J,'dwt'}) doest the same for the DWT (|FWT|)
%   filterbank tree.
%
%   The function internally calls WFBTINIT and passes wt and all
%   additional parameters to it.
%
%   Examples:
%   ---------
%
%   The following two examples create a multirate identity filterbank
%   using a tree of depth 3. In the first example, the filterbank is
%   identical to the DWT tree:
%
%      [g,a] = wfbt2filterbank({'db10',3,'dwt'});
%      filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   In the second example, the filterbank is identical to the full
%   wavelet tree:
%
%      [g,a] = wfbt2filterbank({'db10',3,'full'});
%      filterbankfreqz(g,a,1024,'plot','linabs','posfreq');
%
%   See also: wfbtinit
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/wfbt2filterbank.php

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


complainif_notenoughargs(nargin,1,'WFBT2FILTERBANK');

% build the tree
wt = wfbtinit({'strict',wt},varargin{:});

% Pick just nodes with outputs
wtPath = 1:numel(wt.nodes);
wtPath(noOfNodeOutputs(1:numel(wt.nodes),wt)==0)=[];

rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,wt);


