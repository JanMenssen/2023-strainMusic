function res = strainChunkLSQSE(iter,iAngle,res,grid,additional)
% STRAINCHUNKLSQSE       calculates the strain from displacements using LSQSE
%
%  this function calculates the stran parameters using a linear LSQ
%  estimations.
%  The fields <sxx>, <sxy>, <syy> and <syx> are added to the RES struct
%
%  For 3D this routine is speeded up by using chunks of data, the chunksize
%  is a field in the additional parameter
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
%   22-nov-2018  JM   initial settings, adapted from <strainLSDSE>


%% check input arguments

  narginchk(4,5);
  nargoutchk(1,1);
  
  if nargin == 4, additional= []; end
  assert(isfield(additional,'window'),'StrainMusic:strainLSQSE','wrong additional parameter');
  assert(isfield(additional,'ncores'),'StrainMusic:strainLSQSE','ncores parameter unknown');

  window = additional.window;
  numWorkers = getNumWorkers(additional.ncores);

  % chunksize 
  
  chunksize = [];
  if isfield(additional,'chunksize'), chunksize = additional.chunksize; end
  
  % get real world displacements
  
  [ax,lat,ele] = getRealWorldDisplacement(res{iAngle},iter,grid);
    
%% determine function, depending on 2D/3D

  switch length(grid.size)
  
     
    case 2   
      [szz,szx,sxx,sxz] = strainLSQSE_2D(res{iAngle}.dispindx{iter},ax,lat,grid,window,numWorkers);      
    case 3       
      [szz,szx,szy,sxx,sxz,sxy,syy,syz,syx] = strainLSQSE_3D(res{iAngle}.dispindx{iter},ax,lat,ele,grid,window,numWorkers,chunksize);
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

function  [szz,szx,szy,sxx,sxz,sxy,syy,syz,syx] = strainLSQSE_3D(indices,axdata,latdata,eledata,grid,window,ncores,chunksize)

  % allocate memory for results
  
  sxx = zeros(1,length(indices));
  sxy = zeros(1,length(indices));  
  sxz = zeros(1,length(indices));
  
  syx = zeros(1,length(indices));
  syy = zeros(1,length(indices));  
  syz = zeros(1,length(indices));
    
  szx = zeros(1,length(indices));
  szy = zeros(1,length(indices));  
  szz = zeros(1,length(indices));
 
  % determine the number of chunkds and the starting and ending indices
  
  if isempty(chunksize), chunksize = length(indices); end
  
  nChunks = ceil(length(indices)/chunksize);
  chunkStartIndx = 1:chunksize:(nChunks*chunksize);
  chunkEndIndx = chunkStartIndx + chunksize - 1;
  chunkEndIndx(chunkEndIndx>length(indices)) = length(indices);
  
  % get the positions of the indices to perform the LSQ and sort these data
  % original a transpose was done, this results in an sorting error, now
  % the transpose is done just before lsese estimation

  xPos = grid.x(indices);
  yPos = grid.y(indices);
  zPos = grid.z(indices);

  [xPos,indxSorted] = sort(xPos);
  yPos = yPos(indxSorted);
  zPos = zPos(indxSorted);
  
  % sort the displacementes, originals are not needed anymore, so they are
  % overwritten by sorted values
  
  axdata = axdata(indxSorted);
  latdata = latdata(indxSorted);
  eledata = eledata(indxSorted);
  
  for iChunk=1:nChunks
    
    indxInChunk = chunkStartIndx(iChunk):chunkEndIndx(iChunk);
    [indxPosition,indxInData] = indxIn3Dchunk(indxInChunk,xPos,yPos,zPos,window);
    curChunkSize = length(indxInChunk);
    
    % data to workers
    
    xPosToWorkers = xPos(indxPosition);
    yPosToWorkers = yPos(indxPosition);
    zPosToWorkers = zPos(indxPosition);
    
    axDataToWorkers = axdata(indxPosition);
    latDataToWorkers = latdata(indxPosition);
    eleDataToWorkers = eledata(indxPosition);
    
    % allocate memory
    
    sxxWorkers = NaN(1,curChunkSize);
    sxyWorkers = NaN(1,curChunkSize);
    sxzWorkers = NaN(1,curChunkSize);
    
    syxWorkers = NaN(1,curChunkSize);
    syyWorkers = NaN(1,curChunkSize);
    syzWorkers = NaN(1,curChunkSize);
  
    szxWorkers = NaN(1,curChunkSize);
    szyWorkers = NaN(1,curChunkSize);
    szzWorkers = NaN(1,curChunkSize);
    
    % and now for the current chunk in parallel

    parfor (indx = 1:curChunkSize,ncores)
 
      indxIn = indxInData(indx);
      
      % set the displacement indices that are within the kernel on TRUE
   
      inZkernel = (zPosToWorkers >= zPosToWorkers(indxIn)-window(1) & zPosToWorkers <= zPosToWorkers(indxIn)+window(1));
      inXkernel = (xPosToWorkers >= xPosToWorkers(indxIn)-window(2) & xPosToWorkers <= xPosToWorkers(indxIn)+window(2));
      inYkernel = (yPosToWorkers >= yPosToWorkers(indxIn)-window(3) & yPosToWorkers <= yPosToWorkers(indxIn)+window(3));
      inKernel = inXkernel & inZkernel & inYkernel;
    
      nrPoints = sum(inKernel);
      if (nrPoints > 1)

        % create matrix for coordinates in kernel, transpose is done because
        % it's not doen due to sorting data
      
        axPoints = zPosToWorkers(inKernel)';
        latPoints = xPosToWorkers(inKernel)';
        elePoints = yPosToWorkers(inKernel)';
      
        A = [axPoints latPoints elePoints ones(nrPoints,1)];
        B = (A'*A)\A';

        % in Axial direction

        data = axDataToWorkers(inKernel)';
          
        szzWorkers(indx) = B(1,:) * data;
        szxWorkers(indx) = B(2,:) * data;
        szyWorkers(indx) = B(3,:) * data;
      
        % in lateral direction

        data = latDataToWorkers(inKernel)';
            
        sxzWorkers(indx) = B(1,:) * data;
        sxxWorkers(indx) = B(2,:) * data;
        sxyWorkers(indx) = B(3,:) * data;
      
        % elevational directrion
      
        data = eleDataToWorkers(inKernel)';
      
        syzWorkers(indx) = B(1,:) * data;
        syxWorkers(indx) = B(2,:) * data;
        syyWorkers(indx) = B(3,:) * data;
     
      end
    end
    
    % set values calculated for chunk in whole vector
    
    sxx(indxInChunk) = sxxWorkers;
    sxy(indxInChunk) = sxyWorkers;
    sxz(indxInChunk) = sxzWorkers;

    syx(indxInChunk) = syxWorkers;
    syy(indxInChunk) = syyWorkers;
    syz(indxInChunk) = syzWorkers;

    szx(indxInChunk) = szxWorkers;
    szy(indxInChunk) = szyWorkers;
    szz(indxInChunk) = szzWorkers;

  end
  
  sxx(indxSorted) = sxx;
  sxy(indxSorted) = sxy;
  sxz(indxSorted) = sxz;

  syx(indxSorted) = syx;
  syy(indxSorted) = syy;
  syz(indxSorted) = syz;

  szx(indxSorted) = szx;
  szy(indxSorted) = szy;
  szz(indxSorted) = szz;

end
  

%% indxIn3Dchunk
%
%   function retuns the indices of the data that is needed for one chunk

function [indxPos,indxData] = indxIn3Dchunk(indx,xPos,yPos,zPos,kernel)

  xMax = max(xPos(indx));  xMin = min(xPos(indx));
  yMax = max(yPos(indx));  yMin = min(yPos(indx));
  zMax = max(zPos(indx));  zMin = min(zPos(indx));
  
  inZkernel = (zPos >= zMin-kernel(1) & zPos <= zMax+kernel(1));
  inXkernel = (xPos >= xMin-kernel(2) & xPos <= xMax+kernel(2));
  inYkernel = (yPos >= yMin-kernel(3) & yPos <= yMax+kernel(3));
  
  indxPos = inZkernel & inXkernel & inYkernel;
  
  indxTmp = false(size(indxPos));
  indxTmp(indx) = true;
  indxData = find(indxTmp(indxPos));
  
end
