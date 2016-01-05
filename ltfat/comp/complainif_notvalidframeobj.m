function complainif_notvalidframeobj(F,callfun)

if nargin<2
    callfun = mfilename;
end

if ~isstruct(F) || ~isfield(F,'frana')
  error('%s: Agument F must be a frame definition structure.',...
        upper(callfun));
end;


%
%   Url: http://ltfat.sourceforge.net/doc/comp/complainif_notvalidframeobj.php

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

