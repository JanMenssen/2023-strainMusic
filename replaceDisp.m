function res = replaceDisp(iter,angle,res,~,inStruct)
% REPLACEDISP      replace calculated displacement
%
%   this routine replaces the displacements for iteration <iter>
%   and angle <angle> with the displacment found in the file. This
%   file is an additonal parameter
%
%     syntax : res = replaceDisp(iter,angle,res,grid,inStruct)
%
%   with
%     - iter        : current iteration
%     - angle       : current angle
%     - res         : resulst structure that is modifed
%     - inStruct    : structure with the field <filename>
%
%  see also

%   Modifications
%     30-mar-2016  JM     initial version

%% paramater handling

  narginchk(5,5);
  nargoutchk(1,1);
     
  assert(~isempty(inStruct),'StrainMusic:PostProcessing','no optional parameters'); 
  assert(isfield(inStruct,'file'),'StrainMusic:PostProcessing','file parameter unknown');
  
%% try to open the file and check fields

  s = load(inStruct.file);
  assert(isfield(s.RES{angle},'axf'),'StrainMusic:PostProcessing','wrong RES file');
  assert(isfield(s.RES{angle},'latf'),'StrainMusic:PostProcessing','wrong RES file');
  assert(isfield(s.RES{angle},'elef'),'StrainMusic:PostProcessing','wrong RES file');
  
%% and replace the displacements

  res{angle}.axf{iter} = s.RES{angle}.axf{end};
  res{angle}.latf{iter} = s.RES{angle}.latf{end};
  res{angle}.elef{iter} = s.RES{angle}.elef{end};

%% and roi and displacement indices

  res{angle}.roi{iter} = s.RES{angle}.roi{end};
  res{angle}.dispindx{iter} = s.RES{angle}.dispindx{end};
  
end  
    