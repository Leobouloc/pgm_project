% LTFAT - Filterbanks
%
%  Peter L. Søndergaard, 2011 - 2014
%
%  Transforms and basic routines
%    FILTERBANK             - Filter bank
%    UFILTERBANK            - Uniform Filter bank
%    IFILTERBANK            - Inverse normal/uniform filter bank
%    FILTERBANKWIN          - Evaluate filterbank window
%    FILTERBANKLENGTH       - Length of filterbank to expand signal
%    FILTERBANKLENGTHCOEF   - Length of filterbank to expand coefficients
%
%  Auditory inspired filterbanks
%    CQT                    - Constant Q transform
%    ICQT                   - Inverse constant Q transform
%    ERBLETT                - Erb-let transform
%    IERBLETT               - Inverse Erb-let transform
%    INSDGFB                - Say what?
%
%  Filter generators
%    CQTFILTERS             - Logaritmically spaced filters
%    ERBFILTERS             - ERB-spaced filters
%  
%  Window construction and bounds
%    FILTERBANKDUAL         - Canonical dual filters
%    FILTERBANKTIGHT        - Canonical tight filters
%    FILTERBANKREALDUAL     - Canonical dual filters for real signals
%    FILTERBANKREALTIGHT    - Canonical tight filters for real signals
%    FILTERBANKBOUNDS       - Frame bounds of filter bank
%    FILTERBANKREALBOUNDS   - Frame bounds of filter bank for real signals
%    FILTERBANKRESPONSE     - Total frequency response (a frame property)
%
%  Auxilary
%    FILTERBANKFREQZ        - Frequency responses of filters
%    NONU2UFILTERBANK       - Non-uni. to uniform filterbank transformation
%    U2NONUCFMT             - Change format of coefficients
%    NONU2UCFMT             - Change format of coefficients back
%
%  Plots
%    PLOTFILTERBANK         - Plot normal/uniform filter bank coefficients
%
%  For help, bug reports, suggestions etc. please send email to
%  ltfat-help@lists.sourceforge.net
%
%   Url: http://ltfat.sourceforge.net/doc/filterbank/Contents.php

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


