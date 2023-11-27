function [template,kernel,res] = rigidBlockMatch(iter,angle,res,pre,post,inStruct)
% RIGIDBLOCKMATCH             finds kernel/templates using rigid method
%
%   this function finds all possible kernel/template combination that can
%   be used displacment estimation using normalized cross-correlation.
% 
%     syntax [template,kernel] = rigidBlockMatch(iter,angle,res,pre,post,inputs)
%
%   with
%       iter     : iteration number
%       res      : results structure
%       pre      : dataframe before compression
%       post     : dataframe after compression
%       inStruct : contains additonal fields (windowsize and max
%                : displacement
%
%   see also

% Modifications
%   15-jul-2015   JM    initial version
%   12-aug-2015   JM    now for 2D and 3D
%   13-aug-2015   JM    cell arrays are returned
%   20-aug-2015   JM    bug found when outside data, res now output arg
%   15-feb-2016   JM    dispIndxInUse changed to dispindx{iter}
%   12-jan-2017   JM    other interface due to 161215-01
%   01-aug-2017   JM    pnts and offset now calculated in new calcIndices (170801-01)
%   30-jul-2020   GH    fixed when nrElem=0 (all indices out frame removed)

%% parameter handling

%-jm  narginchk(6,6);
%-jm  nargoutchk(3,3);
     
  assert(~isempty(inStruct),'StrainMusic:BlockMatch','no optional parameters'); 
  assert(isfield(inStruct,'window'),'StrainMusic:BlockMatch','window parameter unknown');
  assert(isfield(inStruct,'max_disp'),'StrainMusic:BlockMatch','window parameter unknown');

%% some small calculations
  
  kernelSize = inStruct.window + inStruct.max_disp;
  templateSize = inStruct.window;

%% get the points and offsets that are used for creating templates kernels

  [nrElem,axPnts,latPnts,elePnts,offset_ax,offset_lat,offset_ele,res.dispindx{iter}] = calcIndices(iter,angle,res,size(pre),kernelSize);
  
%% and now define all combinations  
  kernel.size = 2 * kernelSize + 1;
%   kernel.data{nrElem} = []; 
  kernel.data = cell(1,nrElem);
  
  template.size = 2 * templateSize + 1;
%   template.data{nrElem} = [];
  template.data = cell(1,nrElem);
  
  switch ndims(pre)
    
    case 2
      for i=1:nrElem
        template.data{i} = pre(axPnts(i)-templateSize(1):axPnts(i)+templateSize(1),latPnts(i)-templateSize(2):latPnts(i)+templateSize(2));
        kernel.data{i} = post(axPnts(i)-kernelSize(1)+offset_ax(i):axPnts(i)+kernelSize(1)+offset_ax(i),latPnts(i)-kernelSize(2)+offset_lat(i):latPnts(i)+kernelSize(2)+offset_lat(i));
      end
       
    case 3
      for i=1:nrElem
        template.data{i} = pre(axPnts(i)-templateSize(1):axPnts(i)+templateSize(1),latPnts(i)-templateSize(2):latPnts(i)+templateSize(2),elePnts(i)-templateSize(3):elePnts(i)+templateSize(3));
        kernel.data{i} = post(axPnts(i)-kernelSize(1)+offset_ax(i):axPnts(i)+kernelSize(1)+offset_ax(i),latPnts(i)-kernelSize(2)+offset_lat(i):latPnts(i)+kernelSize(2)+offset_lat(i),elePnts(i)-kernelSize(3)+offset_ele(i):elePnts(i)+kernelSize(3)+offset_ele(i));
      end
      
    otherwise
      error('StrainMusic:BlockMatch','wrong dimensions');
  end

  
end  
  
