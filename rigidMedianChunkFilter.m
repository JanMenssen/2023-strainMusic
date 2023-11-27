function data = rigidMedianChunkFilter(data,indices,grid,inStruct)
% RIGIDMEDIANCHUNKFILTER          performs a median filer
%
%   this is a median filter, usign a rigid kernel and filters only the diaplacment 
%   points that are inside the kernel. Parallel implementation
%   Kernel size must be defined 
%
%   This filter divides the points that should be filtered in a number of
%   chunks, defined by the chunksize. For each chunk, a rigid Median filter
%   is performed parallel. Using thsi method is faster, less data transfer
%     
%     syntax [data,outStruct] = rigidMedianChunkFilter(data,indices,grid,inStruct)
%
%   with
%     - data        : vector of data points to filter
%     - indices     : indices of these data points
%     - grid        : structure as defined in USGRID, with x, y and z
%                   : vector containing coordinates of indices
%     - inStruct    : structure with field <kernel> (optional)
%
%   see also RIGIDMEDIANFILTERS

%   Modifications
%      21-nov-2018  JM    initial version, adapted from RigidMedianFilter

%% argument handling
  
  narginchk(4,4);
  nargoutchk(1,1);
  
  assert(~isempty(inStruct),'StrainMusic:rigidMedianChunkFilter','kernel not given');
  assert(isfield(inStruct,'kernel'),'StrainMusic:rigidMedianChunkFilter','kernel not given'); 
  assert(isfield(inStruct,'ncores'),'StrainMusic:rigidMedianChunkFilter','ncores parameter unknown');

  % check rotation, because we rotate the grid instead of the kernel, the
  % angle becomes negative
  
  rotationAngle = [];
  if isfield(inStruct,'angle'), rotationAngle = -inStruct.angle; end
  
  % check the chunkSize
  
  chunkSize = [];
  if isfield(inStruct,'chunksize'), chunkSize = inStruct.chunksize; end
  
%% handle the number of Workers used

   numWorkers = getNumWorkers(inStruct.ncores);

%% get all coordinates for the indices

  xPosTmp = [];
  yPosTmp = [];
  zPosTmp = [];
 
  if ~isempty(grid.x), xPosTmp = grid.x(indices); end  
  if ~isempty(grid.y), yPosTmp = grid.y(indices); end  
  if ~isempty(grid.z), zPosTmp = grid.z(indices); end
  
%% and rotate the grid if angle is given

  if isempty(rotationAngle)
    xPos = xPosTmp;
    yPos = yPosTmp;
    zPos = zPosTmp;
  else
    xPos = xPosTmp*cosd(rotationAngle) - zPosTmp*sind(rotationAngle);
    yPos = yPosTmp;
    zPos = xPosTmp*sind(rotationAngle) + zPosTmp*cosd(rotationAngle);
  end
  
  clear xPosTmp yPosTmp zPosTmp;
  
%% check 2D/3D and set function for it

  switch length(grid.size)
    
    case 2
       indxfunc = @indxIn2Dkernel;
       indxInChunkfunc = @indxIn2Dchunk;
    
    case 3
       indxfunc = @indxIn3Dkernel;
       indxInChunkfunc = @indxIn3Dchunk;
    
    otherwise
      error('StrainMusic:rigidMedianFilter','wrong number of directions')
  
  end
  
  kernel = inStruct.kernel;
  
%% determine the number of chunks and the starting/ending index for each chunk

  if isempty(chunkSize), chunkSize = length(indices); end
  
  nChunks = ceil(length(indices)/chunkSize);
  chunkStartIndx = 1:chunkSize:(nChunks*chunkSize);
  chunkEndIndx = chunkStartIndx + chunkSize - 1;
  chunkEndIndx(chunkEndIndx>length(indices)) = length(indices);

%% sort to reduce the data send to the kernel, sorting is done in axial direction

  yPosSorted = [];
  zPosSorted = [];
  
  [xPosSorted,indxSorted] = sort(xPos);
  if ~isempty(yPos), yPosSorted = yPos(indxSorted); end
  if ~isempty(zPos), zPosSorted = zPos(indxSorted); end
  
  dataSorted = data(indxSorted);
  
%% now for every chunk

  xPosToWorkers = [];
  yPosToWorkers = [];
  zPosToWorkers = [];
  
  for iChunk=1:nChunks
   
    % calculate the indices for this chunk
    
    indxInChunk = chunkStartIndx(iChunk):chunkEndIndx(iChunk);
    [indxPosition,indxInData] = feval(indxInChunkfunc,indxInChunk,xPosSorted,yPosSorted,zPosSorted,kernel);
    curChunkSize = length(indxInChunk);
    
    % grid and data to filter send to the workers
    
    if ~isempty(xPosSorted), xPosToWorkers = xPosSorted(indxPosition); end
    if ~isempty(yPosSorted), yPosToWorkers = yPosSorted(indxPosition); end
    if ~isempty(zPosSorted), zPosToWorkers = zPosSorted(indxPosition); end
    
    dataToWorker = dataSorted(indxPosition);
    tmpData = NaN(1,curChunkSize);
        
    % and perform in parallel  
   
    parfor (indx=1:curChunkSize,numWorkers)

      indxToWorker = indxInData(indx);
    	inKernel = feval(indxfunc,indxToWorker,xPosToWorkers,yPosToWorkers,zPosToWorkers,kernel);
      tmpData(indx) = nanmedian(dataToWorker(inKernel));  
  
    end
    data(indxInChunk) = tmpData;

  end
  
  % done, reshuffle data to original position
  
  data(indxSorted) = data;
  
end

%% indx2Dkernel
%
%   indices that are within the kernel are true (2D)

function indices = indxIn2Dkernel(indx,xPos,~,zPos,kernel)
   
  inZkernel = (zPos >= zPos(indx)-kernel(1) & zPos <= zPos(indx)+kernel(1));
  inXkernel = (xPos >= xPos(indx)-kernel(2) & xPos <= xPos(indx)+kernel(2));

  indices = inZkernel & inXkernel; 

end


%% Indx3Dkernel
%
%   indices that are within the kernel are true (2D)

function indices = indxIn3Dkernel(indx,xPos,yPos,zPos,kernel)

  inZkernel = (zPos >= zPos(indx)-kernel(1) & zPos <= zPos(indx)+kernel(1));
  inXkernel = (xPos >= xPos(indx)-kernel(2) & xPos <= xPos(indx)+kernel(2));
  inYkernel = (yPos >= yPos(indx)-kernel(3) & yPos <= yPos(indx)+kernel(3));
  
  indices = inZkernel & inXkernel & inYkernel; 

end


%% indxIn2Dchunk
%
%   function retuns the indices of the data that is needed for one chunk

function [indxPos,indxData] = indxIn2Dchunk(indx,xPos,~,zPos,kernel)

  xMax = max(xPos(indx));  xMin = min(xPos(indx));
  zMax = max(zPos(indx));  zMin = min(zPos(indx));
  
  inZkernel = (zPos >= zMin-kernel(1) & zPos <= zMax+kernel(1));
  inXkernel = (xPos >= xMin-kernel(2) & xPos <= xMax+kernel(2));
  
  indxPos = inZkernel & inXkernel;
  
  indxTmp = false(size(indxPos));
  indxTmp(indx) = true;
  indxData = find(indxTmp(indxPos));
    
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

  
