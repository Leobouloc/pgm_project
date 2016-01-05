function [AF,BF]=uwpfbtbounds(wt,L,varargin)
%UWPFBTBOUNDS Frame bounds of Undecimated WPFBT
%   Usage: fcond=uwpfbtbounds(wt,L);
%          [A,B]=uwpfbtbounds(wt,L);
%
%   UWPFBTBOUNDS(wt,L) calculates the ratio B/A of the frame bounds
%   of the undecimated wavelet packet filterbank specified by wt for a 
%   system of length L. The ratio is a measure of the stability of the 
%   system. 
%
%   [A,B]=uwfbtbounds(wt,L) returns the lower and upper frame bounds
%   explicitly. 
%
%   See WFBT for explanation of parameter wt.
%
%   See also: uwpfbt, filterbankbounds
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/uwpfbtbounds.php

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


complainif_notenoughargs(nargin,2,'UWPFBTBOUNDS');

definput.flags.scaling={'sqrt','scale','noscale'};
definput.flags.interscaling = {'intsqrt', 'intscale', 'intnoscale'};
[flags]=ltfatarghelper({},definput,varargin);

wt = wfbtinit({'strict',wt},'nat');

for ii=1:numel(wt.nodes)
   a = wt.nodes{ii}.a;
   assert(all(a==a(1)),sprintf(['%s: One of the basic wavelet ',...
                                'filterbanks is not uniform.'],...
                                upper(mfilename)));
end

% Do the equivalent filterbank using multirate identity property
[gmultid,amultid] = wpfbt2filterbank(wt,flags.interscaling);

% Scale filters 
gmultid = comp_filterbankscale(gmultid, amultid, flags.scaling);

if nargout<2
   AF = filterbankbounds(gmultid,1,L);
elseif nargout == 2
   [AF, BF] = filterbankbounds(gmultid,1,L);
end

