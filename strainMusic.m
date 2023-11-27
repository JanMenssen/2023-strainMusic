function res = strainMusic(preFile,postFile,paramFile,roiFile,gridFile,resFile)
%% STRAINMUSIC     main function for the strain algorithm
%
%   strainMusic is the main function of the strain (=displacement) algorithm 
%   used in the MUSIC group. Displacement is calculated between 2 US frames (pre 
%   and post) using the settings in a parameter file
%
%     input
%       - preFile   : name of file of first US frame that is used
%       - postFile  : name of file of second US frame that is used
%       - paramFile : name of the file with parameters for calculation
%       - roiFile   : name of the file containing the ROI to be used
%       - gridFile  : name of the file containing the used grid 
%       - resFile   : name of the file the results are stored
%     output
%       - res       : standard structure with results of the displacment and strain
%
%   roiFile, gridFile and resFile are optional parameters
%
%   see also

%   Modifications
%     05-jun-2015   JM   initial version
%     05-okt-2015   JM   bug found in prepareData  
%     07-okt-2015   JM   a lot of changes in <calcDisplacment> due to memory
%     15-mar-2016   JM   now works for planes in 3D
%     15-apr-2016   JM   bugs solved if no displacement is calculated
%     09-may-2016   JM   nr of Runs in calcDisplacement now in parameter file
%                        user info changed
%     12-jan-2017   JM   scatterinterpolant only calculated when needed (161215-01)
%     04-may-2017   JM   green code and copyright in header set to 2017
%                        now prepareData is faster, allocation only in enveloppe
%     08-jun-2017   JM   bug 170515-01 solved, variables in memory in <calcDisplacement>
%     11-jun-2018   JM   after if-end, comma deleted, copyright changed
%     21-nov-2018   JM   for all if-end, comma deleted
%     27-jun-2019   JM   bug 190627-01, if only one disp-indices left over (line 335)

%% opening screen

  showInfo;
  
%% check input and output arguments

  narginchk(3,6);
  nargoutchk(0,1);
  
  if (nargin < 6), resFile = NaN; end        % use NaN if no resFile is required
  if (nargin < 5), gridFile = []; end
  if (nargin < 4), roiFile = []; end
  
%% read the pres and post frame, header is used from the 

  disp('reading data and other required files ...');
  
  [preFrameData,frameHdr] = readUSdataFrame(preFile);
  [postFrameData,~] = readUSdataFrame(postFile);

%% read the parameter file and fill fields that are not defined

  [param,nrAngles] = readStrainParam(paramFile);
  
%% if ROI file is given, read it, else create a ROI matrix based on the data

  if (~isempty(roiFile))
    roi = readROIfile(roiFile);
  else
    roi = createROIfromData(preFrameData);
  end
     
%% read the gridFile, if not given use the grid from the header

  if (~isempty(gridFile))
    [usGrid,dispGrid] = readGridFile(gridFile);
  else
    [usGrid,dispGrid] = createGrid(preFrameData,frameHdr);
  end

%% create RES structure

  res = createRESstruct(nrAngles,param);   

%% some preprocessing, iteration independent

  res{1} = setVersion(res{1});
  res{1} = setFileNames(res{1},preFile,postFile,roiFile,gridFile,paramFile);  
  res = initROI(res,param{1}.no_iter,roi,dispGrid);

%% and clear some variables that are not needed anymore

  clearvars roi dispGrid;

%% starting parallel pool if not already started

  if isempty(gcp), parpool; end
  
%% now peform the algorithm for each iteration and each angle
%
%       - calculate displacements
%       - filter the displacements 
%       - do some postproccesing
%
%  the <angle> loop is inside these routines

  for iterNr=1:param{1}.no_iter

    disp(' ');
    disp(['iteration ' num2str(iterNr,'%1d') ' of ' num2str(param{1}.no_iter,'%1d')]);
    
    % calculate the displacement, filter it and calculate strain (postprocessing)

    res = calcDisplacement(iterNr,res,param,preFrameData,postFrameData);
    res = filterDisplacement(iterNr,res,param,usGrid);
    res = postProcessing(iterNr,res,param,usGrid);
      
  end
  
%% done, store the RES structure if required

  for iAngle = 1:nrAngles
    if param{iAngle}.final_iter_only, res{iAngle} = removeResPrevIterations(param{iAngle}.no_iter,res{iAngle}); end
    if param{iAngle}.filtered_disps_only, res{iAngle} = removeUnfilteredResults(param{iAngle}.no_iter,res{iAngle}); end
  end
  
  saveTheResult(resFile,res);

  disp('  ');
  if ~isempty(gcp), delete(gcp); end
  disp('done ...');
  disp('  ');
  
end


%% calcDisplacement
%
%   this function calculates the displacments. Estimateted displacements
%   are stored in the field <ax>, <lat> en <ele> for axial, lateral and
%   elevational displacement

function res = calcDisplacement(iter,res,param,pre,post)
    
  nrChar = 0;
  nrAngles = size(res,2);
     
  for iAngle=1:nrAngles

    % first clear the used variables, previous values were kept in memory
    % which results in a bug
    
    clear cmax kernel template peaks;
    
    % prepare the data, hilbert or raw and make doubles
    
    [preData,postData] = prepareData(iter,pre{iAngle},post{iAngle},param{iAngle});
    [n,res{iAngle}] = adaptNrDispIndices(iter,res{iAngle},'init',param{iAngle}.nCCF{iter});
    
    for iRun = 1:n

      % user info
      
      fprintf(1,repmat('\b',1,nrChar));
      nrChar = fprintf('estimating displacements (if required) ...     [%d/%d - %d/%d]',iAngle,nrAngles,iRun,n);
      
      [~,res{iAngle}] = adaptNrDispIndices(iter,res{iAngle},'begin');
    
      % find the kernels/windows combinations for CCF 

      if ~isempty(param{iAngle}.ccf_blockmatch(iter).func)
        [template,kernel,res{iAngle}] = param{iAngle}.ccf_blockmatch(iter).func(iter,iAngle,res{iAngle},preData,postData,param{iAngle}.ccf_blockmatch(iter).additional);       
      end
      
      % perform cross correlation and find peaks. Note : using this approach disables the functionality to find
      % better peaks using more angles

      if ~isempty(param{iAngle}.ccf_method(iter).func)
        [ccf,cmax,peaks] = param{iAngle}.ccf_method(iter).func(template,kernel,param{iAngle}.ccf_method(iter).additional);
      end
      
      % and do a better peak finding

      if ~isempty(param{iAngle}.ccf_peakinterp(iter).func)
        [cmax,peaks] = param{iAngle}.ccf_peakinterp(iter).func(ccf,peaks,param{iAngle}.ccf_peakinterp(iter).additional);
      end   

      % and store peaks in <res> struct if peaks are found, offset_xxx is empty
      % if that direction doesn't exists

      if exist('cmax','var')
        res{iAngle}.cmax{iter} = [res{iAngle}.cmax{iter} cmax];
      end
      
      if (exist('kernel','var') && exist('template','var'))
        offset = round((kernel.size - template.size + 1)/2);
        [offset_ax,offset_lat,offset_ele] = getOffsets(iter,iAngle,res{iAngle},[]);
      end
      
      if exist('peaks','var') 
        if ~isempty(peaks.ax), res{iAngle}.ax{iter} = [res{iAngle}.ax{iter} (peaks.ax + offset_ax - offset(1))]; end
        if ~isempty(peaks.lat), res{iAngle}.lat{iter} = [res{iAngle}.lat{iter} (peaks.lat + offset_lat - offset(2))]; end
        if ~isempty(peaks.ele), res{iAngle}.ele{iter} = [res{iAngle}.ele{iter} (peaks.ele + offset_ele - offset(3))]; end  
      end
      
      [~,res{iAngle}] = adaptNrDispIndices(iter,res{iAngle},'end');
    end
    
  [~,res{iAngle}] = adaptNrDispIndices(iter,res{iAngle},'done');    
  end

  % and be prepared for next user info message
  
  fprintf('\n');
  
end


%% filterDisplacement
%
%   this function filters the displacements as calculated by <calcDisplament>. 
%   The user specifies which filters are used. filtered values are stored
%   in <axf> <latf> and <elef> fields of the res-struct

function res = filterDisplacement(iter,res,param,usgrid)

  nrAngles = size(res,2);
  
  for iAngle = 1:nrAngles
    
    fprintf('filtering displacements (if required) ...      [%d/%d]\n',iAngle,nrAngles);
  
    indices = res{iAngle}.dispindx{iter};
    grid = usgrid{iAngle};
    
    % axial filter
    
    data = [];
    if ~isempty(res{iAngle}.ax), data = res{iAngle}.ax{iter}; end
    if ~isempty(param{iAngle}.axfilter(iter).func)
      data = param{iAngle}.axfilter(iter).func(data,indices,grid,param{iAngle}.axfilter(iter).additional);
    end
    res{iAngle}.axf{iter} = data; 
    
    % lateral filter
    
    data = [];
    if ~isempty(res{iAngle}.lat), data = res{iAngle}.lat{iter}; end
    if ~isempty(param{iAngle}.latfilter(iter).func)
      data = param{iAngle}.latfilter(iter).func(data,indices,grid,param{iAngle}.latfilter(iter).additional);
    end
    res{iAngle}.latf{iter} = data; 

    % elevational filter
    
    data = [];
    if ~isempty(res{iAngle}.ele), data = res{iAngle}.ele{iter}; end
    if ~isempty(param{iAngle}.elefilter(iter).func)
      data = param{iAngle}.elefilter(iter).func(data,indices,grid,param{iAngle}.elefilter(iter).additional);
    end
    res{iAngle}.elef{iter} = data;
    
  end
  
end


%% postProcessing
%
%     this function performs some post-processing if specified by the user.

function res = postProcessing(iter,res,param,usgrid)
  
  nrAngles = size(res,2);
  for iAngle = 1:nrAngles
  
    grid = usgrid{iAngle};
    nrSteps = size(param{iAngle}.postprocessing,1);
    for iStep=1:nrSteps

      fprintf('postprocessing (if required) ...               [%d/%d - %d/%d]\n',iAngle,nrAngles,iStep,nrSteps);
      
      if ~isempty(param{iAngle}.postprocessing(iStep,iter).func)       
        res = param{iAngle}.postprocessing(iStep,iter).func(iter,iAngle,res,grid,param{iAngle}.postprocessing(iStep,iter).additional);  
      end
      
    end    
  end
end


%% adaptNrDispIndices
%
%   this function splits the number of displacement indices in parts so the
%   displacement calculation can be done in parts so less memory is needed

function [n,res] = adaptNrDispIndices(iter,res,mode,nRun)
  
  persistent firstIndx;
  persistent lastIndx;
  persistent ntmp;
  
%JM  ntalAtaTime  = 250;
    
  switch lower(mode)
  
    % initialisation,  
    
    case 'init'
 
      res.alldispindx = res.dispindx{iter};
      ntal = length(res.alldispindx);
     
      if isempty(nRun), nRun = ntal; end
               
      firstIndx = 1;
%JM      lastIndx = ntalAtaTime;
      lastIndx = nRun;
      if (lastIndx >= ntal), lastIndx = ntal; end
        
%JM      ntmp = ceil(ntal/ntalAtaTime);
      ntmp = ceil(ntal/nRun);
      
    % prepare for the next run    
    
    case 'begin'
      
      res.dispindx{iter} = res.alldispindx(firstIndx:lastIndx);   
      
    % end, adapt variables  
    
    case 'end'

      res.alldispindx = [res.alldispindx(1:firstIndx-1) res.dispindx{iter} res.alldispindx(lastIndx+1:end)];
      firstIndx = firstIndx + length(res.dispindx{iter});
      lastIndx = lastIndx + length(res.dispindx{iter});

      ntal = length(res.alldispindx);
%-jm bug 190726-01    if (firstIndx >= ntal), firstIndx = []; end
      if (firstIndx > ntal), firstIndx = []; end
      if (lastIndx >= ntal), lastIndx = ntal; end
    
    % done, copy dispindx and remove tempory field
    
    case 'done'
       
      res.dispindx{iter} = res.alldispindx;
      res = rmfield(res,'alldispindx');

    otherwise
      error('MATLAB:ADAPNTDISPINDICES','wrong mode');
  
  end    

  % and return the largest number of 
  
  n = ntmp;
  
end


%% prepareData
%
%   ths function prepares the framedata so it can be used in the next steps
%  

function [outPre,outPost] = prepareData(iter,inPre,inPost,param)

  % allocate memory
  
%JM  outPre = NaN * ones(size(inPre),'double');
%JM  outPost = NaN * ones(size(outPre),'double');
  
  switch lower(param.signaltype{iter})
    
    % RAW : convert only to double
    
    case 'raw'
      
      outPre = double(inPre);
      outPost = double(inPost);
     
    % ENV : create enveloppe data
    
    case 'env'

      % allocate memory 
      
      outPre = NaN * ones(size(inPre),'double');
      outPost = NaN * ones(size(outPre),'double');
    
      [~,nrLateral,nrElevational] = size(inPre);
     
      for k=1:nrElevational
        for j=1:nrLateral
          outPre(:,j,k) = abs(hilbert(inPre(:,j,k)));
          outPost(:,j,k) = abs(hilbert(inPost(:,j,k)));
        end
      end
        
    otherwise
      error('StrainMusic:StrainMusic','signaltype wrong');
  end  
     
end


%% removeResPreviousIteration
%
%   this function removes the reulst from previous iterations to save
%   memory

function res = removeResPrevIterations(no_iter,res)
   
  for iter=1:(no_iter-1)
      
    % fields that always exists

    res.ax{iter} = [];
    res.lat{iter} = [];
    res.ele{iter} = [];

    res.axf{iter} = [];
    res.latf{iter} = [];
    res.elef{iter} = [];

    res.cmax{iter} = [];
    
    % fields that are optional
  
    if isfield(res,'sxx'), res.sxx{iter} = []; end
    if isfield(res,'sxy'), res.sxy{iter} = []; end
    if isfield(res,'sxz'), res.sxz{iter} = []; end
    if isfield(res,'syx'), res.syx{iter} = []; end
    if isfield(res,'syy'), res.syy{iter} = []; end
    if isfield(res,'sxz'), res.sxz{iter} = []; end
    if isfield(res,'syz'), res.syz{iter} = []; end
    if isfield(res,'szz'), res.szz{iter} = []; end
 
  end  
end


%% removeUnfilteredResults
%
%   this function replaces the unfiltered results with an empty value

function res = removeUnfilteredResults(no_iter,res)

  for iter=1:no_iter
    res.ax{iter} = [];
    res.lat{iter} = [];
    res.ele{iter} = [];
  end  
end


%% saveTheResult
%
%  save the result structure if required

function saveTheResult(resFile,RES) %#ok<INUSD>
  
  if ~isnan(resFile), save(resFile,'RES'); end

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
    if ~isempty(pathstr), filename = sprintf('%s%s%s',pathstr,filesep,filename); end
    s = load(filename);
    assert(isfield(s,'USHEADER'),'StrainMusic:readUSdataFrame','wrong header file');
    hdr = s.USHEADER;
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


%% readStrainParam
%
%   this routine reads the paramater structure that describes how the displacement 
%   and strain is calculated

function [param,nrAngles] = readStrainParam(paramFile)

  assert(exist(paramFile,'file')==2,'StrainMusic:readStrainParam','File doesn''t exist');
  s= load(paramFile);
  assert(isfield(s,'PARAM'),'StrainMusic:readStrainParam','not a parameter file');
  param = s.PARAM;
  
  nrAngles = size(param,2);
  
end


%% createROIfomData
%
%   this function creates a ROI structure using the entire frame

function roi = createROIfromData(data)

  nrAngles = size(data,2);
  roi{nrAngles} = [];
  
  for iAngle=1:nrAngles, roi{iAngle} = true(size(data{iAngle})); end

end


%% createGrid
%
%     creates a grid, using the defaults (whatever they maybe)

function [usgrid,dispgrid] = createGrid(data,hdr)
 
  nrAngles = size(data,2);
  
  [nrAxialPnts,nrLateralPnts,nrElevationalPnts] = size(data{1});
  switch ndims(data)
    
    case 2
       
      axialDistances = hdr.c/(2*hdr.fs) * (1:nrAxialPnts);
      axPoints = repmat(axialDistances,1,nrLateralPnts,nrElevationalPnts);
      axPoints = axPoints(:);
  
      lateralDistances = hdr.pitch(1) * (1:nrLateralPnts);
      lateralDistances = lateralDistances - mean(lateralDistances);
      latPoints = repmat(lateralDistances,nrAxialPnts,1,nrElevationalPnts);
      latPoints = latPoints(:);
      
      elePoints = [];
     
    case 3
      
      axialDistances = hdr.c/(2*hdr.fs) * (1:nrAxialPnts);
      axPoints = repmat(axialDistances,1,nrLateralPnts,nrElevationalPnts);
      axPoints = axPoints(:);
  
      lateralDistances = hdr.pitch(1) * (1:nrLateralPnts);
      lateralDistances = lateralDistances - mean(lateralDistances);
      latPoints = repmat(lateralDistances,nrAxialPnts,1,nrElevationalPoints);
      latPoints = latPoints(:);

      elevationDistances = hdr.pitch(2) * (1:nrElevationalPnts);
      elevationDistances = elevationDistances - mean(elevationDistances);
      elePoints = repmat(elevationDistances,nrAxialPnts,nrLateralPnts,1);
      elePoints = elePoints(:);
      
    otherwise
      error('StrainMusic:StrainMusic','wrong dimension');
  end
  
  usgrid{nrAngles} = [];
  dispgrid{1,nrAngles} = [];

  for iAngle=1:nrAngles
    
    usgrid{iAngle}.x = (latPoints + sind(hdr.xmitAngles(iAngle))*axPoints)';
    usgrid{iAngle}.y = elePoints;  
    usgrid{iAngle}.z = (axPoints * cosd(hdr.xmitAngles(iAngle)))';
    usgrid{iAngle}.size = size(data);  
  
    dispgrid{1,iAngle}.indices = 1:numel(data{iAngle});
  
  end
  
end  


%% createRESstruct
%
%   creates an emtpy RES structure

function RES = createRESstruct(nrAngles,param)

  RES = cell(1,nrAngles);
  
  for iAngle=1:nrAngles
    
    % the displacment indices
    
    RES{iAngle}.dispindx = [];
    
    % fields containing the ROI values for each displacement pixel per
    % iteration
  
    RES{iAngle}.roi = [];
   
    % displacment fields 
    for iter=1:param{iAngle}.no_iter
    
      % unfiltered displacement for all iteration
  
      RES{iAngle}.ax{iter} = [];
      RES{iAngle}.lat{iter} = [];
      RES{iAngle}.ele{iter} = [];
  
      % filtered displacement results per iteration
  
      RES{iAngle}.axf{iter} = [];
      RES{iAngle}.latf{iter} = [];
      RES{iAngle}.elef{iter} = [];
  
      % fields containing the value fo the peak of the cross-correlion fucntion
      % for each iteration
  
      RES{iAngle}.cmax{iter} = [];
 
    end 
 
    % for every angle a specific parameter file can be used
    
    RES{iAngle}.paramfile = [];
  
  end  
  
  % fields refering to files, data etc to be used in the calculations
  
  RES{1}.prefile = [];
  RES{1}.postfile = [];
  RES{1}.roifile = [];
  RES{1}.gridfile = [];
  
  % field referint to version of strain software, OS and Matlab version
  % (strings)
  
  RES{1}.strain_version = [];
  RES{1}.os_version = [];
  RES{1}.os_system = [];
  RES{1}.matlab_version = [];
    
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


%% setFileNames
%
%   function sets the names of the used file in the RES structure

function res = setFileNames(res,prefile,postfile,roifile,gridfile,paramfile)

  res.prefile = prefile;
  res.postfile = postfile;
  res.roifile = roifile;
  res.gridfile = gridfile;
  res.paramfile = paramfile;
  
end

%% setVersion
%
%   this function fills the version fields in the RES structure with
%   information

function res = setVersion(res)

  res.strain_version = 2.0;
  res.os_system = computer;
  
  if ispc, [~,res.os_version] = system('ver'); end
  
  if ismac, [~,res.os_version] = system('sw_vers'); end
  
  if isunix
  end
  
  res.matlab_version = version;

end


%% showInfo
%
%   function shows program information

function showInfo()

disp('                                                                           ');
disp('    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
disp('    *                                                                     *');
disp('    *                       S t r a i n M u s i c                         *');
disp('    *                                                                     *');
disp('    *            a 2D/3D strain algorithm for ultrasound data             *');
disp('    *                                                                     *');
disp('    *     Rik Hansen & Jan Menssen - Medical UltraSound Imaging Center    *');
disp('    *       Radboud University Nijmegen Medical Center - 2015-2019(c)     *');
disp('    *                                                                     *');
disp('    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *');
disp('                                                                           ');

end
