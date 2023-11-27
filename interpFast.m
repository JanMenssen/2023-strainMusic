function [ccfMax,peaks] = interpFast (ccfmat,peaks,inStruct)
%INTERPFAST       peforms a 1/2/3D fit for better peak finding
%
%   Given the CCF matrix and the peaks found in this matrix, this routine
%   does a more exact peak finding using 1-, 2- or 3-D peak interpolation
%       
%       syntax: [ccfMax,peaks] = interpFast(ccfmat,peaks,inStruct)
%
%   with
%      - ccfmat.data{iPeak} : cross-correlation matrix (for all frames)
%      - peaks              : peaks for all frames
%      - inStruct           : struct with field <additional> containing 
%                             <fieldsnIter>, <psize>, <method>, <ncores>
%                           <nIter> number of iterative steps to 
%                               interpolate on finer grid around the peak.
%                           <psize> number of ccf points around 
%                               peak(1sided) in x,y,z for peak interpolation. 
%                               Set to [] to use all points
%                           <method> interpolation method: 'nearest',
%                               'linear',spline', 'cubic'
%                           <ncores> number of cores for parallel processes.
%
%     outputs
%         ccfMax  : maximum value of ccf
%         peaks   : peaks found after peak interpolation
%
% NOTE: if psize <2: problems can occur if using spline (too few data
% points if one corner of ccfmat).
%
% Framework addapted from J. Menssen (parpar_2D) and S. Fekkes (iterative
% intrp. algorithm).

%  Modfications:
%       01-dec-2015 GH  Original version
%       29-jan-2016 GH  if peak.ax/lat/ele isempty -> ax/lat/eletop = 1 .
%       04-mar-2016 GH  NaNs removed in the ccfparts (error if peak is in
%                       corner ccfmat).
%       05-mar-2016 GH  Bug fix at border: no extrapolation but zeros at 
%                       points outside ccf (this is border of max disp)
%       12-apr-2016 GH  Update ccfmat to ccfmat.data for new alg. version
%       08-jun-2016 GH  Update in comments

%% argument handling

  narginchk(3,3);
  nargoutchk(2,2);
  
  assert(~isempty(inStruct),'StrainMusic:interpFast','additional parameters not given');
  assert(isfield(inStruct,'nIter'),'StrainMusic:interpFast','nIter not given'); 
  assert(isfield(inStruct,'method'),'StrainMusic:interpFast','method not given'); 
  assert(isfield(inStruct,'ncores'),'StrainMusic:interpFast','ncores parameter unknown');
  assert(isfield(inStruct,'psize'),'StrainMusic:interpFast','psize parameter unknown');
  
%% handle the number of Workers used

  numWorkers = getNumWorkers(inStruct.ncores);
   
%% check 1D, 2D or 3D interpolation function
 
  switch numel(find(ccfmat.size>1))
    case 1
      intpfunc = @intp1D;
    case 2
      intpfunc = @intp2D;    
    case 3
      intpfunc = @intp3D;
    otherwise
      error('StrainMusic:interpFast','wrong number of dimensions')
  end
  
  
%% Define current top and max before interpolation
  
  nrPeaks = length(ccfmat.data);
  ccfMax = NaN * ones(1,nrPeaks);
  
  if isempty(peaks.ax) 
    axTop(1,1:nrPeaks) = 1;
  else
    axTop = peaks.ax(1,:);
  end
  
  if isempty(peaks.lat) 
    latTop(1,1:nrPeaks) = 1;
  else
    latTop = peaks.lat(1,:);
  end
  
  if isempty(peaks.ele) 
    eleTop(1,1:nrPeaks) = 1;
  else
    eleTop = peaks.ele(1,:);
  end

%% initialiation of interpolation function

  % create parts of ccfmat.data around peak (squared, psize, one-sided)
  
  psize = inStruct.psize;
  ccfParts = feval(@createCCFparts, ccfmat.data, axTop, latTop, eleTop, psize);
 
  % settings intpfunc:
  nIter = inStruct.nIter;
  method = inStruct.method;
  
%% perform interpolation
     
  parfor (iPeak = 1:nrPeaks,numWorkers)
    [axShift, latShift, eleShift, ccfMax(iPeak)] = feval(intpfunc,ccfParts{iPeak},nIter,method);
    axTop(iPeak) = axTop(iPeak) + axShift;
    latTop(iPeak) = latTop(iPeak) + latShift;
    eleTop(iPeak) = eleTop(iPeak) + eleShift;
  end

  % Set new peaks after interpolation
  
  if ~isempty(peaks.ax);  peaks.ax(:) = axTop;   end
  if ~isempty(peaks.lat); peaks.lat(:) = latTop; end 
  if ~isempty(peaks.ele); peaks.ele(:) = eleTop; end

  
end

%% createCCFparts
%
% create ccf parts around the peak used for peak interpolation (psize is
% one-sided)

function  ccfParts = createCCFparts(ccfmat, axTop, latTop, eleTop, psize)

  nrPeaks = length(ccfmat); 
  ccfParts = repmat({[]},1,nrPeaks); 

  if ~isempty(psize)
      
    % Create parts around the peak of psize (if possible)
    
    for iPeak = 1:nrPeaks

      axStart = (axTop(iPeak)-psize); if axStart<1; axStart = 1;  end
      axEnd =   (axTop(iPeak)+psize); if axEnd>size(ccfmat{iPeak},1); axEnd = size(ccfmat{iPeak},1); end

      latStart = (latTop(iPeak)-psize); if latStart<1; latStart = 1; end
      latEnd =   (latTop(iPeak)+psize); if latEnd>size(ccfmat{iPeak},2); latEnd = size(ccfmat{iPeak},2); end

      eleStart = (eleTop(iPeak)-psize); if eleStart<1; eleStart = 1; end
      eleEnd =   (eleTop(iPeak)+psize); if eleEnd>size(ccfmat{iPeak},3); eleEnd = size(ccfmat{iPeak},3); end

      ccfParts{iPeak} = ccfmat{iPeak}(axStart:axEnd, latStart:latEnd, eleStart:eleEnd);

      % NaN in ccfmat set to 0 (bad estimate)
      ccfParts{iPeak}(isnan(ccfParts{iPeak}(:))) = 0; 
    end
               
  else
    
    % or use entire ccfmat for peak interpolation
    
    for iPeak = 1:nrPeaks
        ccfParts{iPeak} = ccfmat{iPeak};        
        ccfParts{iPeak}(isnan(ccfParts{iPeak}(:))) = 0;   % NaN in ccfmat set to 0 (bad estimate)
    end
    
  end
end


%% Intp1D
%
%   Iterative 1D interpolation: every iteration peak is interpolated grid
%   around peak, new peak found, new finer grid is generated around new peak
%   ... till final step: location of peak is subsample displacment.

function [rowShift, colShift, elShift, ccfMax] = intp1D(ccf, iter, method)

  [rowMax,colMax, elMax]  = ind2sub(size(ccf),find(ccf==max(ccf(:))));
  [iym]           = (-1:1)';

  rowMaxOrg = rowMax(1);
  colMaxOrg = colMax(1);
  elMaxOrg = elMax(1);

  rowMax = rowMax(1);
  colMax = colMax(1);
  elMax  = elMax(1);

  for i = 1:iter

    iy = iym.*.5^i + rowMax;

    ccfi = interp1(ccf,iy,method,0);

    [rowM,~,~]  = ind2sub(size(ccfi),find(ccfi==max(ccfi(:))));

    rowM = rowM(1);
    rowMax = rowMax + (rowM-2)*.5^i;

  end

  ccfMax   = max(ccfi(:));

  rowShift = rowMax(1)-rowMaxOrg;
  colShift = colMax(1)-colMaxOrg;
  elShift  = elMax (1)-elMaxOrg;

end


%% Intp2D
%
%   Iterative 2D interpolation: every iteration peak is interpolated grid
%   around peak, new peak found, new finer grid is generated around new peak
%   ... till final step: location of peak is subsample displacment.

function [rowShift, colShift, elShift, ccfMax] = intp2D(ccf, iter, method)

  [rowMax,colMax, elMax]  = ind2sub(size(ccf),find(ccf==max(ccf(:))));
  [ixm,iym]           = meshgrid(-1:1:1 ,-1:1:1);

  rowMaxOrg = rowMax(1);
  colMaxOrg = colMax(1);
  elMaxOrg = elMax(1);

  rowMax = rowMax(1);
  colMax = colMax(1);
  elMax  = elMax(1);

  for i = 1:iter
    
    iy = iym.*.5^i + rowMax;
    ix = ixm.*.5^i + colMax;
    
    ccfi = interp2(ccf,ix,iy,method,0); %if outside ccf means; outside max disp so zero!  No interp
    
    [rowM,colM, ~]  = ind2sub(size(ccfi),find(ccfi==max(ccfi(:))));
    
    rowM = rowM(1);
    colM = colM(1);
    
    rowMax = rowMax + (rowM-2)*.5^i;
    colMax = colMax + (colM-2)*.5^i;
       
  end

  ccfMax   = max(ccfi(:));

  rowShift = rowMax(1)-rowMaxOrg;
  colShift = colMax(1)-colMaxOrg;
  elShift  = elMax (1)-elMaxOrg;

end


%% Intp3D
%
% Iterative 3D interpolation: every iteration peak is interpolated grid
% around peak, new peak found, new finer grid is generated around new peak
% ... till final step: location of peak is subsample displacment.

function [rowShift, colShift, elShift, ccfMax] = intp3D(ccf, iter, method)

  [rowMax,colMax, elMax] = ind2sub(size(ccf),find(ccf==max(ccf(:))));
  [ixm,iym,izm] = meshgrid(-1:1:1 ,-1:1:1, -1:1:1);

  rowMaxOrg = rowMax(1);
  colMaxOrg = colMax(1);
  elMaxOrg = elMax(1);

  rowMax = rowMax(1);
  colMax = colMax(1);
  elMax  = elMax(1);
        
  for i = 1:iter

    iy = iym.*.5^i + rowMax;
    ix = ixm.*.5^i + colMax;
    iz = izm.*.5^i + elMax;

    ccfi = interp3(ccf,ix,iy,iz,method,0);

    [rowM,colM, elM]  = ind2sub(size(ccfi),find(ccfi==max(ccfi(:))));

    rowM = rowM(1);
    colM = colM(1);
    elM  = elM(1);

    rowMax = rowMax + (rowM-2)*.5^i;
    colMax = colMax + (colM-2)*.5^i;
    elMax  = elMax  + (elM-2)*.5^i;

  end

  ccfMax   = max(ccfi(:));

  rowShift = rowMax(1)-rowMaxOrg; %Shift by interpolation!
  colShift = colMax(1)-colMaxOrg; %Shift by interpolation!
  elShift  = elMax (1)-elMaxOrg; %Shift by interpolation!

end




  