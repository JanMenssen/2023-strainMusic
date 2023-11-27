function resRaw = interpResult(resEnv,resRaw,usGridEnv,usGridRaw,roiFile,downSample)
% interpDisplacement    interpolates displacements form eveloppe to raw
%
%   this routine interpolates results obtained from the <strainMusic> algorithme 
%   to other x/y/z points. This routine can be used to interpolate results obtained 
%   from the downsampled enveloppe data to raw RF data
%
%     syntax : interpRes(resEnv,usGridEnv,usGridEnv,usGridRaw)
%
%       with
%         - resEnv        ; name of input res-file (enveloppe)
%         - resRaw        : name of res file containing interpolated data
%         - usGridEnv     : name of Grid file used for result (enveloppe)
%         - usGridRaw     : grid file used to calculate interpolated data
%         - roiFIle       : optional ROI file
%         - downSample    : factor used for downsampling and to mulitply
%                         : displacements
%
%   NOTE : currently ONLY the displacement fields are interpolated and ONLY
%        : 3D is supported 

%   Modifications
%      22-mar-2016  JM    initial version
%      08-may-2016  JM    downSample factor added
%      06-jun-2016  JM    bug 160608_01 solved (grid converted to doubles)

  narginchk(4,6);
  nargoutchk(0,1);
  if (nargin <= 4), roiFile = []; end;
  if (nargin <= 5), downSample = []; end;
  if isempty(downSample), downSample = [1 1 1]; end;
    
%% read the files and display indices for the data to be interpolated

  [envGrid,~] = readGridFile(usGridEnv);
  [rawGrid,rawDisp] = readGridFile(usGridRaw);
  
  if (~isempty(roiFile))
    roi = readROIfile(roiFile);
  else
    roi{size(rawGrid,2)} = [];
    for iAngle = 1:size(rawGrid,2), roi{iAngle} = true(rawGrid{iAngle}.size); end;
  end
  
  oldres = load(resEnv);
  oldres = oldres.RES;

%% get indices for last and next iteration

  nrAngles = size(oldres,2);
  nrIter = size(oldres{1}.dispindx,2);
  
  RES{nrAngles} = [];
  RES = initROI(RES,1,roi,rawDisp);
  
%% and interpolate

  for iAngle = 1:nrAngles
    
    envIndices = oldres{iAngle}.dispindx{nrIter};
    rawIndices = RES{iAngle}.dispindx{1};
  
    % displacements
    
    data = [];
    if ~isempty(oldres{iAngle}.axf{nrIter})
      data = oldres{iAngle}.axf{nrIter} * downSample(1);
      data = interpolator3D(data,envIndices,envGrid{iAngle},rawIndices,rawGrid{iAngle});
    end
    RES{iAngle}.axf{1} = data;
    
    data = [];
    if ~isempty(oldres{iAngle}.latf{nrIter})
      data = oldres{iAngle}.latf{nrIter} * downSample(2);
      data = interpolator3D(data,envIndices,envGrid{iAngle},rawIndices,rawGrid{iAngle});
    end
    RES{iAngle}.latf{1} = data;

    data = [];
    if ~isempty(oldres{iAngle}.elef{nrIter})
      data = oldres{iAngle}.elef{nrIter} * downSample(3);
      data = interpolator3D(data,envIndices,envGrid{iAngle},rawIndices,rawGrid{iAngle});
    end
    RES{iAngle}.elef{1} = data;
    
  end
  
%% and store result file

  save(resRaw,'RES');
  
end


%% interpolator
%
%       interpolates data

function outdata = interpolator3D(indata,inIndices,inGrid,outIndices,outGrid)

  % get all x,y and z indices and check they are a in a plane
  
  prevXpoints = double(inGrid.x(inIndices))';
  prevYpoints = double(inGrid.y(inIndices))';
  prevZpoints = double(inGrid.z(inIndices))';
  
  newXpoints = double(outGrid.x(outIndices))';
  newYpoints = double(outGrid.y(outIndices))';
  newZpoints = double(outGrid.z(outIndices))';
  
  xEqual = isequal(min(prevXpoints),max(prevXpoints));
  yEqual = isequal(min(prevYpoints),max(prevYpoints));
  zEqual = isequal(min(prevZpoints),max(prevZpoints));
  
  if (xEqual || yEqual || zEqual)
  
    if xEqual
      f = scatteredInterpolant(prevYpoints,prevZpoints,indata','linear','nearest');
      outdata = f(newYpoints,newZpoints);      
    end
     
    if yEqual
      f = scatteredInterpolant(prevXpoints,prevZpoints,indata','linear','nearest');
      outdata = f(newXpoints,newZpoints);
    end
    
    if zEqual
      f = scatteredInterpolant(prevXpoints,prevYpoints,indata','linear','nearest');
      outdata = f(newXpoints,newYpoints);
    end
    
  else
    
    f = scatteredInterpolant(prevXpoints,prevYpoints,prevZpoints,indata','linear','nearest');
    outdata = f(newXpoints,newYpoints,newZpoints);
    
  end
  
end


%% readGridFile
%
%   this function read the GRID file and returns both fields USGRID en
%   DISPGRID

function [usGrid,dispGrid] = readGridFile(filename)

  s = load(filename);
  assert(isfield(s,'USGRID'),'StrainMusic:readGridFile','not a GRID file');
  assert(isfield(s,'DISPGRID'),'StrainMusic:readGridFIle','not a GRID file');
  
  usGrid = s.USGRID;
  dispGrid = s.DISPGRID;
  
end


%% initROI
%
%  This routine finds the indices in the ROI which are true and combines
%  these with the indices given in the dispgrid
 
function res = initROI(res,nrIterations,roi,dspGrd)

  nrAngles = size(res,2);
  
  iterInGrid = size(dspGrd,1);
  angleInGrid = size(dspGrd,2);

  for iAngle=1:nrAngles
       
    if angleInGrid == 1
      dspAngle = 1;
    else  
      dspAngle = iAngle;
    end
    
    roiIndices = find(roi{iAngle});

    for iter = 1:nrIterations
      
      res{1,iAngle}.roi{iter}.size = size(roi{iAngle});
          
      if iterInGrid == 1
        dspIter = 1;
      else 
        dspIter = iter;
      end
      
      % only displacmente indices inside ROI
      
      res{iAngle}.dispindx{iter} = uint32(intersect(dspGrd{dspIter,dspAngle}.indices,roiIndices'));

    end   
  end  
end


%% readROIfile
%
%   this routome reads the ROI file given by filename. The ROI is stored "old-syle"
%   using a matrix or "new-style", using a structure. In both cases the "new-style"
%   structure is returned

function roi = readROIfile(filename)

  narginchk(1,1);
  nargoutchk(1,1);
  
  s = load(filename);
  assert(isfield(s,'ROI'),'StrainMusic:readROIfile','not a ROI file');
  if ~iscell(s.ROI)
    s.ROI(isnan(s.ROI)) = 0;
    roi{1} = logical(s.ROI);
  else
    roi = s.ROI;
  end
  
end

