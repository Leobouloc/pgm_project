function coef=framenative2coef(F,coef);
%FRAMENATIVE2COEF  Convert coefficient from native format
%   Usage: cout=framenative2coef(F,cin);
%
%   FRAMENATIVE2COEF(F,coef) convert frame coefficients from the native
%   format of the transform into the common column format.
%
%   See also: frame, framecoef2native
%
%   Url: http://ltfat.sourceforge.net/doc/frames/framenative2coef.php

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
  
complainif_notenoughargs(nargin,2,'FRAMENATIVE2COEF');
complainif_notvalidframeobj(F,'FRAMENATIVE2COEF');

% .native2coef is not a mandatory field
if isfield(F,'native2coef')
   coef=F.native2coef(coef);
end
