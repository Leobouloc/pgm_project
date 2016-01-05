function [AF,BF]=uwfbtbounds(wt,L,varargin)
%UWFBTBOUNDS Frame bounds of Undecimated WFBT
%   Usage: fcond=uwfbtbounds(wt,L);
%          [A,B]=uwfbtbounds(wt,L);
%
%   UWFBTBOUNDS(wt,L) calculates the ratio B/A of the frame bounds
%   of the undecimated filterbank specified by wt for a system of length
%   L. The ratio is a measure of the stability of the system.
%
%   wfbtbounds({w,J,'dwt'},L) calculates the ratio B/A of the frame
%   bounds of the DWT (|FWT|) filterbank specified by w and J for a
%   system of length L.
%
%   [A,B]=UWFBTBOUNDS(...,L) returns the lower and upper frame bounds
%   explicitly.
%
%   See WFBT for explanation of parameter wt and FWT for explanation
%   of parameters w and J.
%
%   See also: uwfbt, filterbankbounds
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/uwfbtbounds.php

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


complainif_notenoughargs(nargin,2,'UWFBTBOUNDS');

definput.flags.scaling={'sqrt','scale','noscale'};
[flags]=ltfatarghelper({},definput,varargin);

wt = wfbtinit({'strict',wt},'nat');

for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),sprintf(['%s: One of the basic wavelet ',...
                                'filterbanks is not uniform.'],...
                                upper(mfilename)));
end

% Do the equivalent filterbank using multirate identity property
[gmultid, amultid] = wfbt2filterbank(wt);

% Scale filters
gmultid = comp_filterbankscale(gmultid, amultid, flags.scaling);

if nargout<2
   AF = filterbankbounds(gmultid,1,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gmultid,1,L);
end

