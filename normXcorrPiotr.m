function [ccf,maxVal,maxPos] = normXcorrPiotr(template,kernel,additional)
%NORMXCORRPIOTRTEST     normalized cross correlation using Piotr Dollar's toolbox
%
%  this function calculatest the displacements between 2 frames with a 
%  normalixed cross-correlation function. This is a parallel implementation
%
%      syntax : [ccf,maxVal,maxPos] = normXcorrPiotr(kernel,template,additional)
%
%  input parameters
%         - kernel     : 3D or 4D matrix containing all kernels. last
%                      : dimension contains the number of kernels
%         - template   : 3/4D matrix containing all templates, last dimension 
%                      : contains the number of templates
%         - additional : additional parameters (.ncores)
%
%  output parameters
%         - ccf       : crosscorrelaton matrix 3D/4D, last dimension is
%                     : number of template/kernel combination
%         - maxValue  : maximum value
%         - maxPos    : location of the peaks

%   Modifications
%        24-jun-2015    JM   initial version
%        28-jul-2015    JM   changed after input/output modification
%        13-aug-2015    JM   normxcorrn used (Piotr Dollar's toolbox)ndm
%        07-okt-2015    JM   now row vectors
%        14-mar-2016    JM   adapted for planes
%        30-jul-2020    GH   fixed when nrElem=0 (all indices out frame removed)

%% some parameter checking

  narginchk(3,3);
  nargoutchk(2,3);

  nrDimensions = size(kernel.size,2);
  nrFrames = size(kernel.data,2);

  ccfsize = kernel.size - template.size + 1;
  
%% addtional parameters ok ?

  assert(isfield(additional,'ncores'),'StrainMusic:NormXcorrPiotr','ncores parameter unknown');

%% handle the number of Workers used

  numWorkers = getNumWorkers(additional.ncores);

%% allocation of some data
%   ccfdata{nrFrames} = [];
  ccfdata = cell(1,nrFrames);

  maxPos.ax = [];
  maxPos.lat = [];
  maxPos.ele = [];
  
  maxVal = NaN * ones(1,nrFrames);
  posMax = NaN * ones(1,nrFrames);

  % do a check if there are plane are if a rotate should be done
  %  - planes : same dimension in kernel and template = 1
  %  - rotate : last dimension in template = 1 (no ele) but kernel is volume
  
  bothPlanes =  find(template.size == 1) == find(kernel.size == 1);  
  shouldRotate = (numel(template.size) == 3) && (template.size(numel(template.size)) == 1) && (kernel.size(numel(kernel.size)) ~= 1);
  
  if shouldRotate

    % normalized cross correlation, however rotate cube/plane
 
    parfor (i=1:nrFrames,numWorkers)
      
      template_tmp = permute(template.data{i},[1 3 2]);
      kernel_tmp = permute(kernel.data{i},[1 3 2]);
      ccf_tmp = normxcorrn(template_tmp,kernel_tmp,'valid');
      ccfdata{i} = ipermute(ccf_tmp,[1 3 2]);
      
      [maxVal(i),tmpPosition] = max(ccfdata{i}(:));
      posMax(i) = tmpPosition(1);    
      
    end
  
  elseif bothPlanes
    
    % normalized crosscorelation for planes, use squeeze and reshape
    
    parfor (i=1:nrFrames,numWorkers)
       
      ccfdata{i} = normxcorrn(squeeze(template.data{i}),squeeze(kernel.data{i}),'valid');
      ccfdata{i} = reshape(ccfdata{i},ccfsize);

      [maxVal(i),tmpPosition] = max(ccfdata{i}(:));
      posMax(i) = tmpPosition(1);
    
    end
  
  else
  
    % normalized crosscorrelation for volumes 
    
    parfor (i=1:nrFrames,numWorkers)
 
      ccfdata{i} = normxcorrn(template.data{i},kernel.data{i},'valid');
 
      [maxVal(i),tmpPosition] = max(ccfdata{i}(:));
      posMax(i) = tmpPosition(1);
    
    end
    
  end
  
%% and convert the index of the top to an ax (1), lateral(2) and elevation (3) position 
  
  switch nrDimensions
    case 2
      [maxPos.ax,maxPos.lat] = ind2sub(ccfsize,posMax);
    case 3
      [maxPos.ax,maxPos.lat,maxPos.ele] = ind2sub(ccfsize,posMax);
    otherwise
      error('StrainMusic:normXcorrPiotr','wrong number of dimensions');
  end

%% and cross correlation matrix 

  ccf.size = ccfsize;
  ccf.data = ccfdata;
  
end
