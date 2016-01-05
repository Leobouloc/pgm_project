function [gu,au,p]=nonu2ufilterbank(g,a)
%NONU2UFILTERBANK   Non-uniform to uniform filterbank transform
%   Usage:  [gu,au]=nonu2ufilterbank(g,a)
%
%   Input parameters:
%         g     : Filters as a cell array of structs.
%         a     : Subsampling factors.
%
%   Output parameters:
%         gu    : Filters as a cell array of structs.
%         au    : Uniform subsampling factor.
%         pk    : Numbers of copies of each filter.
%
%   [gu,au]=NONU2UFILTERBANK(g,a) calculates uniform filterbank gu, 
%   au=lcm(a) which is identical to the (possibly non-uniform) filterbank
%   g,*a in terms of the equal output coefficients. Each filter g{k} 
%   is replaced by p(k)=au/a(k) delayed versions of itself such that
%   z^{-ma(k)}G_k(z) for m=0,...,p-1.
%
%   This allows using the factorisation algorithm when determining
%   filterbank frame bounds in FILTERBANKBOUNDS and
%   FILTERBANKREALBOUNDS and in the computation of the dual filterbank 
%   in FILTERBANKDUAL and FILTERBAKREALDUAL which do not work 
%   with non-uniform filterbanks.
%
%   One can chenge between the coefficient formats of gu, au and 
%   g, a using NONU2UCFMT and U2NONUCFMT in the reverse direction.
%
%   See also: ufilterbank, filterbank, filterbankbounds, filterbankdual
%
%   References:
%     S. Akkarakaran and P. Vaidyanathan. Nonuniform filter banks: New
%     results and open problems. In P. M. C.K. Chui and L. Wuytack, editors,
%     Studies in Computational Mathematics: Beyond Wavelets, volume 10, pages
%     259 -301. Elsevier B.V., 2003.
%     
%
%   Url: http://ltfat.sourceforge.net/doc/filterbank/nonu2ufilterbank.php

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

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

if ~iscell(g) || any(cellfun(@isempty,g)) || ...
   ~all(cellfun(@(gEl) isnumeric(gEl) || ...
   (isstruct(gEl) && (isfield(gEl,'h')||isfield(gEl,'H'))),g))
  error(['%s: g must be a cell array of structs containing filter',...
         ' definition or numeric arrays.'],upper(mfilename));
end;

% Sanitize a
a = comp_filterbank_a(a,numel(g));

if size(a,2)==2 && ~all(a(:,2)==1) && rem(a(:,1),1)~=0
   error(['%s: Filterbanks with fractional subsampling are not',...
          ' supported.'],upper(mfilename)); 
end

% This does lcm(a)
au=filterbanklength(1,a);

% Numbers of copies of each filter
p = au./a(:,1);

if all(a(:,1)==a(1,1))
    % Filterbank is already uniform, there is nothing to be done.
    gu = g;
    au = a;
    return;
end

% Sanitize filters which come as numeric vertors.
gIdx = cellfun(@isnumeric,g);
if any(gIdx)
    g(gIdx) = filterbankwin(g(gIdx),a,'normal');
end

gu=cell(sum(p),1);

auIdx = 1;
for m=1:numel(g)
   for ii=0:p(m)-1
      gu{auIdx} = g{m};
      if(isfield(gu{auIdx},'H'))
         if(~isfield(gu{auIdx},'delay'))
            gu{auIdx}.delay = 0;
         end
         gu{auIdx}.delay = gu{auIdx}.delay-a(m)*ii;
      end
      
      if(isfield(gu{auIdx},'h'))
         if(~isfield(gu{auIdx},'offset'))
            gu{auIdx}.offset = 0;
         end
         gu{auIdx}.offset = gu{auIdx}.offset-a(m)*ii;
      end
      auIdx = auIdx+1;
   end
end



