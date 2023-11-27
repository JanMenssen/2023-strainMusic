function [ccf,maxVal,maxPos] = normXcorrEatonCPU(template,kernel, ~)
%NORMXCORREATONCPU     normalized cross correlation using normxcorrCPU (Daniel Eaton) 
%
%  this function calculatest the displacements between 2 frames with a 
%  normalixed cross-correlation function. Implementation of Eaton's normxcorr, converted to C by Jan Menssen
%
%      syntax : [ccf,maxVal,maxPos] = normXcorrEatonCPU(template,kernel,additional)
%
%  input parameters
%         - kernel     : structure containing 2 fields
%                           - size: size of the kernels. [ax,lat]
%                           - data: cell array of size <number of
%                           kernel-template combinations>, with all kernels
%         - template   : structure containing 2 fields
%                           - size: size of the template. [ax,lat]
%                           - data: cell array of size <number of
%                           kernel-template combinations>, with all templates  
%         - (~/)additional:  none defined for this function
%
%  output parameters
%         - ccf       : structure containing 2 fields
%                           - size: size of the ccf. [ax,lat]
%                           - data: cell array of size <number of
%                           kernel-template combinations>, with all ccfs   
%         - maxValue  : maximum value
%         - maxPos    : location of the maximum value (ie peak) (ax, lat, ele)

%   Modifications
%        08-06-2021   initial version (Anne Saris)

%% some parameter checking

  narginchk(2,3);
  nargoutchk(2,3);

  nrDimensions = size(kernel.size,2);
  nrFrames = size(kernel.data,2);

  ccfsize = kernel.size - template.size + 1;
  
  % check dimensions of kernels (and templates), should be 2D for this function
  assert(nrDimensions == 2, 'wrong number of dimensions')
  
%% prepare input for normXcorrCPU

  seg1 = zeros([template.size nrFrames], 'single');
  seg2 = zeros([kernel.size nrFrames], 'single');
  for i = 1: nrFrames
      seg1(:,:,i) = template.data{i};
      seg2(:,:,i) = kernel.data{i};
  end
  
%% normalized cross-correlation - only for 2D data!

  [~,~,~,ccfAll] = normXcorrCPU(0,seg1(:,:,:),seg2(:,:,:),0); 
  
%% create correct output for strainMUSIC

  posMax = zeros([1, nrFrames], 'single');
  maxVal = zeros([1, nrFrames], 'single');

  ccfdata = cell(1,nrFrames);
  for i = 1: nrFrames
      ccfdata{i} = ccfAll(:,:,i);
      
      [maxVal(i),tmpPosition] = max(ccfdata{i}(:));
      posMax(i) = tmpPosition(1);
  end
 
%% allocation of some data

  maxPos.ax = [];
  maxPos.lat = [];
  maxPos.ele = []; % only required for functionality of strainMUSIC algorithm
    
%% and convert the index of the top to an ax (1), lateral(2) position 
  
  [maxPos.ax,maxPos.lat] = ind2sub(ccfsize,posMax);

%% and cross correlation matrix 

  ccf.size = ccfsize;
  ccf.data = ccfdata;
  
end
