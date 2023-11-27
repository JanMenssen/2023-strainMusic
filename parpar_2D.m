function [ccfMax,peaks] = parpar_2D(ccfmat,peaks,~)
%PARPAR_2D        peforms a parabolic fit for better peak finding
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

  narginchk(2,3);
  nargoutchk(2,3);
  
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
  