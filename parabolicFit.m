  function [ccfMax,peaks] = parabolicFit(ccfmat,peaks,~)
%INTERPPARABOLIC       performs a parabolic fit for better peak finding
%
%   Given the CCF matrix and the peaks found in this matrix, this routine
%   does a more exact peak finding using parabolic fit in all direction
%
%     input
%         ccfmat : cross-correlation matrix (for all frames)
%         peaks  : peaks for all frames
%         inputs : addition, not used input
%
%     outputs
%         ccfMax  : maximumv alue of ccf
%         peaks   : peaks found after parpar
%         outputs : additonal struct as output

%  Modfications
%     01-jul-2015  JM     initial version
%     08-jul-2015  JM     bug found ccfmat.data{iPeak} is 3D (third dim = index)
%                         inputs not used, so removed
%     28-jul-2015  JM     changed after input/output modifications
%     13-aug-2015  JM     changed, because peaks chanbed
%     07-okt-2015  JM     now row-vectors
%     22-mar-2016  JM     adapted for change in ccfmat.data
%     08-jun-2016  JM     now 3D supported and renamed to <interpParabolic>

  narginchk(2,3);
  nargoutchk(2,3);
     
  switch numel(ccfmat.size)

    case 2
      fitFunction = @parabolicFit_2D;
    case 3
      fitFunction = @parabolicFit_3D;
    otherwise
      error('StrainMusic:parabolicFit','wrong number of dimensions')
  end

  [ccfMax,peaks] = fitFunction(ccfmat,peaks);

end


%% parablicFit_2D
%
%   performs a parabolic for 2D data

function [ccfMax,peaks] = parabolicFit_2D(ccfmat,peaks)

  nrPeaks = length(ccfmat.data);
  ccfMax = NaN * ones(1,nrPeaks);
    
  for iPeak = 1:nrPeaks    
    
    xtop = peaks.ax(iPeak);
    ytop = peaks.lat(iPeak);
    
    ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop);
    if ccfMax(iPeak) < 0.9999

      if xtop > 1 && xtop < ccfmat.size(1) && ytop > 1 && ytop < ccfmat.size(2)

        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop))/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop-1,ytop)+ccfmat.data{iPeak}(xtop+1,ytop)) + xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop+1))/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop,ytop-1)+ccfmat.data{iPeak}(xtop,ytop+1)) + ytop;

        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop-1,ytop)+ccfmat.data{iPeak}(xtop+1,ytop))*(ccfmat.data{iPeak}(xtop-1,ytop) -  ...
                         ccfmat.data{iPeak}(xtop+1,ytop))^2 -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop,ytop-1)+ccfmat.data{iPeak}(xtop,ytop+1))*(ccfmat.data{iPeak}(xtop,ytop-1) - ...
                         ccfmat.data{iPeak}(xtop,ytop+1))^2 + ccfmat.data{iPeak}(xtop,ytop);

      elseif xtop > 1 && xtop < ccfmat.size(1)
        
        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop))/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop-1,ytop)+ccfmat.data{iPeak}(xtop+1,ytop)) + xtop;
        peaks.lat(iPeak) = ytop;

        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop-1,ytop)+ccfmat.data{iPeak}(xtop+1,ytop))*(ccfmat.data{iPeak}(xtop-1,ytop)-ccfmat.data{iPeak}(xtop+1,ytop))^2 + ccfmat.data{iPeak}(xtop,ytop);

      elseif ytop > 1 && ytop < ccfmat.size(2)

        peaks.ax(iPeak) = xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop+1))/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop,ytop-1)+ccfmat.data{iPeak}(xtop,ytop+1)) + ytop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop)+ccfmat.data{iPeak}(xtop,ytop-1)+ccfmat.data{iPeak}(xtop,ytop+1))*(ccfmat.data{iPeak}(xtop,ytop-1)-ccfmat.data{iPeak}(xtop,ytop+1))^2 + ccfmat.data{iPeak}(xtop,ytop);

      else
        
        peaks.ax(iPeak)  = xtop;
        peaks.lat(iPeak) = ytop;
        ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop);
        
      end
    end    
  end
   
end

%% parablicFit_3D
%
%   performs a parabolic for 3D data

function [ccfMax,peaks] = parabolicFit_3D(ccfmat,peaks)

  nrPeaks = length(ccfmat.data);
  ccfMax = NaN * ones(1,nrPeaks);

  for iPeak = 1:nrPeaks    
    
    xtop = peaks.ax(iPeak);
    ytop = peaks.lat(iPeak);
    ztop = peaks.ele(iPeak);
    
    ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop,ztop);
    if ccfMax(iPeak) < 0.9999

      if xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && ztop > 1 && ztop < size(ccfmat.data{iPeak},3)

        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
        
      elseif xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && size(ccfmat.data{iPeak},3) < 3 && size(ccfmat.data{iPeak},1) >= 3 && size(ccfmat.data{iPeak},2) >= 3
 
        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
       
      elseif xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && size(ccfmat.data{iPeak},2) < 3 && ztop > 1 && ztop < size(ccfmat.data{iPeak},3) && size(ccfmat.data{iPeak},1) >= 3 && size(ccfmat.data{iPeak},3) >= 3

        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);     
                    
      elseif size(ccfmat.data{iPeak},1) < 3 && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && ztop > 1 && ztop < size(ccfmat.data{iPeak},3) && size(ccfmat.data{iPeak},2) >= 3 && size(ccfmat.data{iPeak},3) >= 3
       
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                    
      elseif size(ccfmat.data{iPeak},1) >= 3 && xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && size(ccfmat.data{iPeak},2) < 3 && size(ccfmat.data{iPeak},3) < 3
       
        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                    
      elseif size(ccfmat.data{iPeak},1) < 3 && size(ccfmat.data{iPeak},2) >= 3 && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && size(ccfmat.data{iPeak},3) < 3
       
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                                     
      elseif size(ccfmat.data{iPeak},1) < 3 && size(ccfmat.data{iPeak},2) < 3 && size(ccfmat.data{iPeak},3) >= 3 && ztop > 1 && ztop < size(ccfmat.data{iPeak},3)
        
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                      
      else

        peaks.ax(iPeak)  = xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = ztop;
        ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop,ztop);
        
      end
    end    
  end
  
end

 
  