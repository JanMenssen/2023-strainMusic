function data = rigidMedianFilter(data,indices,grid,inStruct)
% RIGIDMEDIANFILTER          performs a median filer
%
%   this is a median filter, usign a rigid kernel and filters only the diaplacment
%   points that are inside the kernel. Parallel implementation
%   Kernel size must be defined
%
%     syntax [data,outStruct] = rigidMedianFilter(data,indices,grid,inStruct)
%
%   with
%     - data        : vector of data points to filter
%     - indices     : indices of these data points
%     - grid        : structure as defined in USGRID, with x, y and z
%                   : vector containing coordinates of indices
%     - inStruct    : structure with field <kernel> (optional)
%
%   see also

%   Modifications
%      10-jul-2015  JM    initial version
%      05-aug-2015  JM    changed kernel find (faster)
%      07-aug-2015  JM    now proper directions used
%      12-aug-2105  JM    adapted for 3D
%      27-aug-2015  JM    parallel implementation
%      08-jun-2016  JM    now rotation added (algorihm developed by GH)
%      12-dec-2016  GH    bug solved
%      08-nov-2017  JM    median replaced by nanmedian, )see 171024-01)
%      21-nov-2018  JM    if-end, comma deleted
%      19-apr-2022  RH    nanmedian replaced by median(...,'omitnan')

%% argument handling

narginchk(4,4);
nargoutchk(1,1);

assert(~isempty(inStruct),'StrainMusic:rigidMedianFilter','kernel not given');
assert(isfield(inStruct,'kernel'),'StrainMusic:rigidMedianFilter','kernel not given');
assert(isfield(inStruct,'ncores'),'StrainMusic:rigidMedianFilter','ncores parameter unknown');

% check rotation, because we rotate the grid instead of the kernel, the
% angle becomes negative

rotationAngle = [];
if isfield(inStruct,'angle'), rotationAngle = -inStruct.angle; end

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
        
    case 3
        indxfunc = @indxIn3Dkernel;
        
    otherwise
        error('StrainMusic:rigidMedianFilter','wrong number of directions')
        
end

%% and now do for every point

nPoints = length(indices);
tmpData = data;
kernel = inStruct.kernel;

% NOTE kernel = Instruct.kernel and then argument in indxfunc returns an error

parfor (indx=1:nPoints,numWorkers)
    
    % sets the displacement indices that are within the kernel to true and filter
    
    inKernel = feval(indxfunc,indx,xPos,yPos,zPos,kernel); %#ok<FVAL>
    data(indx) = median(tmpData(inKernel),'omitnan');   %#ok<PFBNS>
    
end

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



