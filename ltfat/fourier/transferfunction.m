function H=transferfunction(g,L)
%TRANSFERFUNCTION  The transferfunction of a filter
%   Usage:  H=transferfunction(g,L);
%
%   TRANSFERFUNCTION(g,L) computes the transferfunction of length L*
%   of the filter defined by g.
%
%   See also: pfilt
%
%   Url: http://ltfat.sourceforge.net/doc/fourier/transferfunction.php

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

error(nargchk(2,2,nargin));

[g,info] = comp_fourierwindow(g,L,upper(mfilename));

H=comp_transferfunction(g,L);

