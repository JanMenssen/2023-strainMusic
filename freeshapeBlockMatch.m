function [template,kernel,res] = freeshapeBlockMatch(iter,angle,res,pre,post,inStruct)
% FREESHAPEBLOCKMATCH             finds kernel/template using the freeshape method
%
%   this function finds a;; possible kernel/template combinations thar are
%   used ito calculate displacemente estimation using the normalized
%   cross-correlation algorithm.
%
%     syntax [template,kernel,outStruct] = freeshapeBlockMatch(iter,angle,res,pre,post,inStruct)
%
%   with
%       iter          : iteration number
%       res           : results structure
%       pre           : dataframe before compression
%       post          : dataframe after compression
%       inStruct      : contains additonal fields (windowsize and max
%                     : displacement
%
%   see aldo RIGIDBLOCKMATCH

%  Modification
%   16-jul-2015    JM   initial version
%   12-aug-2015    JM   for 2D adn 3D, cell-array's returned
%   21-aug-2015    JM   now res-struct returned (indices may change)
%   07-okt-2015    JM   changed due to memory problems in 3d
%   15-feb-2016    JM   dispIndxInUse changed to dispindx{iter}
%   12-jan-2107    JM   interface changed due 161215-01

%% parameter handling

  narginchk(6,6);
  nargoutchk(3,3);

  assert(~isempty(inStruct),'StrainMusic:BlockMatch','no optional parameters'); 
  assert(isfield(inStruct,'window'),'StrainMusic:BlockMatch','window parameter unknown');
  assert(isfield(inStruct,'max_disp'),'StrainMusic:BlockMatch','window parameter unknown');
  
%% some small calculations

  kernelSize = inStruct.window + inStruct.max_disp;
  templateSize = inStruct.window;

%% get the points and offsets that are used for creating templates kernels, free shape, whole matrix is needed

  [nrElem,axPnts,latPnts,elePnts,~,~,~,res.dispindx{iter}] = calcIndices(iter,angle,res,size(pre),kernelSize);
  [offset_ax,offset_lat,offset_ele] = getOffsets(iter,angle,res,kernelSize);
  
%% and now define all combinations,

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
        for lat = (latPnts(i)-kernelSize(2)):(latPnts(i)+kernelSize(2))
          axStrt = axPnts(i)-kernelSize(1)+offset_ax(axPnts(i),lat);
          axEnd = axPnts(i)+kernelSize(1)+offset_ax(axPnts(i),lat);
          kernel.data{i}(:,(lat+kernelSize(2)-latPnts(i)+1)) = post(axStrt:axEnd,lat+offset_lat(axPnts(i),latPnts(i)));             
        end
      end 
  
    case 3
       for i=1:nrElem
        template.data{i} = pre(axPnts(i)-templateSize(1):axPnts(i)+templateSize(1),latPnts(i)-templateSize(2):latPnts(i)+templateSize(2),elePnts(i)-templateSize(3):elePnts+templateSize(3));
        for lat = (latPnts(i)-kernelSize(2)):(latPnts(i)+kernelSize(2))
          for ele = (elePnts(i)-kernelSize(3)):(elePnts(i)+kernelSize(3))
            axStrt = axPnts(i)-kernelSize(1)+offset_ax(axPnts(i),lat,ele);
            axEnd = axPnts(i)+kernelSize(1)+offset_ax(axPnts(i),lat,ele);
            kernel.data{i}(:,(lat+kernelSize(2)-latPnts(i)+1),ele+kernelSize(3)-elePnts(i)+1) = post(axStrt:axEnd,lat+offset_lat(axPnts(i),latPnts(i),elePnts(i)),ele+offset_ele(axPnts(i),latPnts(i),elePnts(i))); 
          end  
        end
      end 
      
    otherwise
      error('StrainMusic:BlockMatch','wrong dimensions');
  end
end  


