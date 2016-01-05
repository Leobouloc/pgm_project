function f=iuwfbt(c,par,varargin)
%IUWFBT   Inverse Undecimated Wavelet Filterbank Tree
%   Usage: f = iuwfbt(c,info)
%          f = iuwfbt(c,wt)
%
%   Input parameters:
%         c       : Coefficients stored in L xM matrix.
%         info/wt : Transform parameters struct/Wavelet tree definition.
%
%   Output parameters:
%         f     : Reconstructed data.
%
%   f = IUWFBT(c,info) reconstructs signal f from the coefficients c*
%   using parameters from info struct. both returned by the UWFBT function.
%
%   f = IUWFBT(c,wt) reconstructs signal f from the wavelet coefficients
%   c using the undecimated wavelet filterbank tree described by wt.
%
%   Please see help for UWFBT description of possible formats of wt.
%
%   Filter scaling:
%   ---------------
%
%   As in UWFBT, the function recognizes three flags controlling scaling of 
%   the filters:
%
%      'sqrt'
%               Each filter is scaled by 1/sqrt(a), there a is the hop
%               factor associated with it. If the original filterbank is
%               orthonormal, the overall undecimated transform is a tight
%               frame.
%               This is the default.
%
%      'noscale'
%               Uses filters without scaling.
%
%      'scale'
%               Each filter is scaled by 1/a.
%
%   If 'noscale' is used, 'scale' must have been used in UWFBT (and vice
%   versa) in order to obtain a perfect reconstruction.
%
%   Examples:
%   ---------
%
%   A simple example showing perfect reconstruction using the "full
%   decomposition" wavelet tree:
%
%     f = greasy;
%     J = 6;
%     c = uwfbt(f,{'db8',J,'full'});
%     fhat = iuwfbt(c,{'db8',J,'full'});
%     % The following should give (almost) zero
%     norm(f-fhat)
%
%   See also:  uwfbt, plotwavelets
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/iuwfbt.php

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

complainif_notenoughargs(nargin,2,'IUWFBT');

if(isstruct(par)&&isfield(par,'fname'))
   complainif_toomanyargs(nargin,2,'IUWFBT');

   wt = wfbtinit({'dual',par.wt},par.fOrder);
   scaling = par.scaling;

   % Use the "oposite" scaling
   if strcmp(scaling,'scale')
       scaling = 'noscale';
   elseif strcmp(scaling,'noscale')
       scaling = 'scale';
   end
else
   %% PARSE INPUT
   definput.import = {'wfbtcommon'};
   definput.flags.scaling={'sqrt','scale','noscale'};
   flags=ltfatarghelper({},definput,varargin);

   % Initialize the wavelet tree structure
   wt = wfbtinit(par,flags.forder);
   scaling = flags.scaling;
end


%% ----- step 2 : Prepare input parameters
wtPath = nodesBForder(wt,'rev');
nodesUps = nodeFiltUps(wtPath,wt);
rangeLoc = rangeInLocalOutputs(wtPath,wt);
rangeOut = rangeInOutputs(wtPath,wt);
%% ----- step 3 : Run computation
f = comp_iuwfbt(c,wt.nodes(wtPath),nodesUps,rangeLoc,rangeOut,scaling);

