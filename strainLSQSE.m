function res = strainLSQSE(iter,iAngle,res,grid,additional)
% STRAINLSQSE       calculates the strain from displacements using LSQSE
%
%  this function calculates the stran parameters using a linear LSQ
%  estimations.
%  The fields <sxx>, <sxy>, <syy> and <syx> are added to the RES struct
%
%     syntax : res = strainLSQSE(iter,res,grid,additional)
%
%  with
%   - iter        : iteration 
%   - res         : result structure (input/output)
%   - grid        : usgrid
%   - additional  : structure containing window size
%
% See for the LSQSE method used Kallel et al 1997

% Modifications
%   28-jul-2015  JM   initial settings
%   24-agu-2015  JM   bug in 2D lateral calculations
%   03-sep-2015  JM   in postprocessing info us used for all angles
%   07-okt-2015  JM   now row-vectors (data)
%   08-okt-2015  KG   bug in sxx direction
%   23-nov-2018  JM   comment if-end removed
%   23-nov-2018  JM   bug found in getRealWorldDisplay
%   23-nov-2018  JM   bug found in Strain calculation

%% check input arguments

  narginchk(4,5);
  nargoutchk(1,1);
  
  if nargin == 4, additional= []; end
  assert(isfield(additional,'window'),'StrainMusic:strainLSQSE','wrong additional parameter');
  assert(isfield(additional,'ncores'),'StrainMusic:strainLSQSE','ncores parameter unknown');

  window = additional.window;
  numWorkers = getNumWorkers(additional.ncores);

  % get real world displacements
  
 [ax,lat,ele] = getRealWorldDisplacement(res{iAngle},iter,grid);
    
%% determine function, depending on 2D/3D

  switch length(grid.size)
  
     
    case 2   
      [szz,szx,sxx,sxz] = strainLSQSE_2D(res{iAngle}.dispindx{iter},ax,lat,grid,window,numWorkers);      
    case 3       
      [szz,szx,szy,sxx,sxz,sxy,syy,syz,syx] = strainLSQSE_3D(res{iAngle}.dispindx{iter},ax,lat,ele,grid,window,numWorkers);
    otherwise
      error('StrainMusic:strainLSQSE','wrong number of dimensions');
  end

%% copy to result structure

  if exist('sxx','var'), res{iAngle}.sxx{iter} = sxx; end
  if exist('sxy','var'), res{iAngle}.sxy{iter} = sxy; end
  if exist('sxz','var'), res{iAngle}.sxz{iter} = sxz; end
  if exist('syx','var'), res{iAngle}.syx{iter} = syx; end
  if exist('syy','var'), res{iAngle}.syy{iter} = syy; end
  if exist('syz','var'), res{iAngle}.syz{iter} = syz; end
  if exist('szx','var'), res{iAngle}.szx{iter} = szx; end
  if exist('szy','var'), res{iAngle}.szy{iter} = szy; end
  if exist('szz','var'), res{iAngle}.szz{iter} = szz; end
  
end


%% getRealWorldDisplacment
%
%   function returns the filtered displacments in real world values instead
%   of sample values

function [ax,lat,ele] = getRealWorldDisplacement(res,iter,grid)

  axSpacing = [];
  latSpacing = [];
  eleSpacing = [];
  
  % axial sapcing
   
  diffZ = diff(grid.z);
  diffZ(diffZ < eps) =[];
  if ~isempty(diffZ), axSpacing = diffZ(1); end

  % lateral spacing
   
  diffX = diff(grid.x);
  diffX(diffX < eps) = [];
  if ~isempty(diffX), latSpacing = diffX(1); end
 
  % elevational apacing
  
  diffY = diff(grid.y);
  diffY(diffY < eps) = [];
  if ~isempty(diffY), eleSpacing = diffY(1); end
  
  % and return data
  
  ax = res.axf{iter} * axSpacing;
  lat = res.latf{iter} * latSpacing;
  ele = res.elef{iter} * eleSpacing;

end


%% strain_LSSQSE_2D
%
%   2D version of calculation strain

function [szz,szx,sxx,sxz] = strainLSQSE_2D(indices,axdata,latdata,grid,window,ncores)

  % get the indices that are used and the data

  xPos = grid.x(indices)';
  zPos = grid.z(indices)';
 
  % pre-allocate some array's

  sizeIndices = size(indices);
  
  szz = NaN * zeros(sizeIndices);
  szx = NaN * zeros(sizeIndices);
  sxx = NaN * zeros(sizeIndices);
  sxz = NaN * zeros(sizeIndices);
  
  % and now for each displacement point

  parfor (indx = 1:length(indices),ncores)

    % set the displacement indices that are within the kernel on TRUE
   
    inZkernel = (zPos >= zPos(indx)-window(1) & zPos <= zPos(indx)+window(1)); %#ok<*PFBNS>
    inXkernel = (xPos >= xPos(indx)-window(2) & xPos <= xPos(indx)+window(2));
    inKernel = inXkernel & inZkernel;
    
    nrPoints = sum(inKernel);

    if (nrPoints > 1)

      % create matrix for coordinates in kernel
      
      axPoints = zPos(inKernel);
      latPoints = xPos(inKernel);
      
      A = [axPoints latPoints ones(nrPoints,1)];
      B = (A'*A)\A';
        
      data = axdata(inKernel)';
      
      szz(indx) = B(1,:) * data;
      szx(indx) = B(2,:) * data;
   
      % in lateral direction

      data = latdata(inKernel)';
      
      sxx(indx) = B(2,:) * data;
      sxz(indx) = B(1,:) * data;
       
    end
  end
  
  % and add strain fields to the <res> struct

end


%% strain_LSSQSE_3D
%
%   3D version of calculation strain

function  [szz,szx,szy,sxx,sxz,sxy,syy,syz,syx] = strainLSQSE_3D(indices,axdata,latdata,eledata,grid,window,ncores)

  % get the indices that are used and the data

  xPos = grid.x(indices)';
  yPos = grid.y(indices)';
  zPos = grid.z(indices)';

  % pre-allocate some array's

  sizeIndices = size(indices);

  sxx = NaN * zeros(sizeIndices);
  sxy = NaN * zeros(sizeIndices);
  sxz = NaN * zeros(sizeIndices);
  syy = NaN * zeros(sizeIndices);
  syx = NaN * zeros(sizeIndices);
  syz = NaN * zeros(sizeIndices);
  szz = NaN * zeros(sizeIndices);
  szx = NaN * zeros(sizeIndices);
  szy = NaN * zeros(sizeIndices);
  
  % and now for each displacement point

  parfor (indx = 1:length(indices),ncores)
 
    % set the displacement indices that are within the kernel on TRUE
   
    inZkernel = (zPos >= zPos(indx)-window(1) & zPos <= zPos(indx)+window(1));
    inXkernel = (xPos >= xPos(indx)-window(2) & xPos <= xPos(indx)+window(2));
    inYkernel = (yPos >= yPos(indx)-window(3) & yPos <= yPos(indx)+window(3));
    inKernel = inXkernel & inZkernel & inYkernel;
    
    nrPoints = sum(inKernel);
    if (nrPoints > 1)

      % create matrix for coordinates in kernel
      
      axPoints = zPos(inKernel);
      latPoints = xPos(inKernel);
      elePoints = yPos(inKernel);
      
      A = [axPoints latPoints elePoints ones(nrPoints,1)];
      B = (A'*A)\A';

      % in Axial direction

      data = axdata(inKernel)';
          
      szz(indx) = B(1,:) * data;
      szx(indx) = B(2,:) * data;
      szy(indx) = B(3,:) * data;
      
      % in lateral direction

      data = latdata(inKernel)';
            
      sxz(indx) = B(1,:) * data;
      sxx(indx) = B(2,:) * data;
      sxy(indx) = B(3,:) * data;
      
      % elevational directrion
      
      data = eledata(inKernel)';
      
      syz(indx) = B(1,:) * data;
      syx(indx) = B(2,:) * data;
      syy(indx) = B(3,:) * data;
     
    end
  end
  
end