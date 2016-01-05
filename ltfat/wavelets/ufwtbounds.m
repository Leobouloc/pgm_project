function [AF,BF]=ufwtbounds(w,J,L,varargin)
%UFWTBOUNDS Frame bounds of Undecimated DWT
%   Usage: fcond=ufwtbounds(w,J,L);
%          [A,B]=ufwtbounds(w,J,L);
%
%   UFWTBOUNDS(w,a,L) calculates the ratio B/A of the frame bounds
%   of the filterbank specified by w and J for a system of length
%   L. The ratio is a measure of the stability of the system.
%
%   [A,B]=UFWTBOUNDS(w,J,L) returns the lower and upper frame bounds
%   explicitly.
%
%   See UFWT for explanation of parameters w and J.
%
%   See also: ufwt, filterbankbounds
%
%   Url: http://ltfat.sourceforge.net/doc/wavelets/ufwtbounds.php

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


complainif_notenoughargs(nargin,3,'UFWTBOUNDS');

if nargout<2
   AF = uwfbtbounds({w,J,'dwt'},L,varargin{:});
elseif nargout == 2
   [AF, BF] = uwfbtbounds({w,J,'dwt'},L,varargin{:});
end

