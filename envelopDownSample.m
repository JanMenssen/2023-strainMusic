function envelopDownSample(dataFilePre,dataFilePost,gridFile,roiFile,factor)
% ENVELOPDOWNSAMPLE
%
%   downsample routine used in the strain software, all necessary
%   files are downsampled by a factor. This factor can be different
%   in axial, lateral and elevational direction
%
%     intput 
%       - datFilePre    : datafile before deformationa
%       - dataFilePost  : datafile after deformation
%       - gridFile      : grid file
%       - roiFile       : roi file
%       - factor        : downsample factor in 3d eg. [5 2 1]
%
%  the name of the downsampled files are the same as the original
%  file postfixed with _downsampled

% Modifications
%     30-mar-2016 JM   initial version

  narginchk(5,5);
  nargoutchk(0,0);

%% read the data files, create an enveloppe file and downsample it

  [USDATA,~] = readUSdataFrame(dataFilePre);
  nrAngles = size(USDATA,2);
  for iAngle = 1:nrAngles
    
    [~,nrLateral,nrElevational] = size(USDATA{iAngle});
     
    for k=1:nrElevational
      for j=1:nrLateral
        USDATA{iAngle}(:,j,k) = abs(hilbert(USDATA{iAngle}(:,j,k)));
      end
    end

    USDATA{iAngle} = USDATA{iAngle}(1:factor(1):end,1:factor(2):end,1:factor(3):end);    
  end  
  save(makeNewName(dataFilePre,1),'USDATA');
  
  [USDATA,~] = readUSdataFrame(dataFilePost);
  nrAngles = size(USDATA,2);
  for iAngle = 1:nrAngles
  
    [~,nrLateral,nrElevational] = size(USDATA{iAngle});
     
    for k=1:nrElevational
      for j=1:nrLateral
        USDATA{iAngle}(:,j,k) = abs(hilbert(USDATA{iAngle}(:,j,k)));
      end
    end

    USDATA{iAngle} = USDATA{iAngle}(1:factor(1):end,1:factor(2):end,1:factor(3):end);    
  end  
  save(makeNewName(dataFilePost,1),'USDATA');

  
%% copy the header file 

  [pathstr,namestr,ext] = fileparts(dataFilePost);
  hdrFileName = sprintf('USHEADER%s%s',namestr(3:(end-6)),ext);
  if ~isempty(pathstr), hdrFileName = sprintf('%s%s%s',pathstr,filesep,hdrFileName); end;

  s = load(hdrFileName);
  assert(isfield(s,'USHEADER'),'StrainMusic:readUSdataFrame','wrong header file');
  USHEADER = s.USHEADER;
  save(makeNewName(hdrFileName,0),'USHEADER');

%% read the ROI file, downsample it

  if ~isempty(roiFile)
    ROI = readROIfile(roiFile);
    nrAngles = size(ROI,2);
    for iAngle = 1:nrAngles
      ROI{iAngle} = ROI{iAngle}(1:factor(1):end,1:factor(2):end,1:factor(3):end);
    end
    save(makeNewName(roiFile,1),'ROI');
  end
  
%% read the grid file and downsample it

  if ~isempty(gridFile)
    
    [USGRID,DISPGRID] = readGridFile(gridFile);
    nrAngles = size(USGRID,2);
    for iAngle = 1:nrAngles

      tmpMatrix = reshape(USGRID{iAngle}.x,USGRID{iAngle}.size);
      tmpMatrix = tmpMatrix(1:factor(1):end,1:factor(2):end,1:factor(3):end);
      USGRID{iAngle}.x = reshape(tmpMatrix,1,numel(tmpMatrix));

      tmpMatrix = reshape(USGRID{iAngle}.y,USGRID{iAngle}.size);
      tmpMatrix = tmpMatrix(1:factor(1):end,1:factor(2):end,1:factor(3):end);
      USGRID{iAngle}.y = reshape(tmpMatrix,1,numel(tmpMatrix));

      tmpMatrix = reshape(USGRID{iAngle}.z,USGRID{iAngle}.size);
      tmpMatrix = tmpMatrix(1:factor(1):end,1:factor(2):end,1:factor(3):end);
      USGRID{iAngle}.z = reshape(tmpMatrix,1,numel(tmpMatrix));
      
      origSize = USGRID{iAngle}.size;
      USGRID{iAngle}.size = size(tmpMatrix);
     
      [axPos,latPos,elePos] = ind2sub(origSize,DISPGRID{1,iAngle}.indices);
      axPos = ceil(axPos/factor(1));
      latPos = ceil(latPos/factor(2));
      elePos = ceil(elePos/factor(3));
      DISPGRID{1,iAngle}.indices = unique(sub2ind(size(tmpMatrix),axPos,latPos,elePos));

    end
    save(makeNewName(gridFile,0),'USGRID','DISPGRID');
  end
end


%% MakeNewName
%
%   creates a filename, same as original but postfixed with _env. For data
%   files (mode = 1) _env is placed before the framenr

function newName = makeNewName(filename,mode)

  switch mode 
  
    case 0
      [path,name,ext] = fileparts(filename);
      newName = [path filesep name '_env' ext];
    
    case 1
      [path,name,ext] = fileparts(filename);
      newName = [path filesep name(1:end-6) '_env' name(end-5:end) ext];
      
    otherwise
      error('StrainMusic:envelopDownSample','wrong mode');
  end
  
end


%% readUSdataFrame
%
%     this routine reads US data. This data is stored in a file, using the 
%     the MUSIC filenname convention. These convention is
%
%        US_<nameofregistration>_<framenr>.mat (5 digits for the framenr)
%
%     Returned is the data matrix and, if required (2nd argument)
%     information regarding the file, stored in a USHEADER struct

function [data,hdr] = readUSdataFrame(filename)

  % load the data file, returned structure should contain field USDATA
  
  s = load(filename);
  assert(isfield(s,'USDATA'),'StrainMusic:readUSdataFrame','wrong file format');
  data = s.USDATA;
  
  % read header file if required
  
  if (nargout == 2) 
    [pathstr,namestr,ext] = fileparts(filename);
    filename = sprintf('USHEADER%s%s',namestr(3:(end-6)),ext);
    if ~isempty(pathstr), filename = sprintf('%s%s%s',pathstr,filesep,filename); end;
    s = load(filename);
    assert(isfield(s,'USHEADER'),'StrainMusic:readUSdataFrame','wrong header file');
    hdr = s.USHEADER;
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