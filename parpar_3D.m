function [ccfMax,peaks] = parpar_3D(ccfmat,peaks,~)
%PARPAR_3D       peforms a parabolic fit for better peak finding
%
%   Given the CCF matrix and the peaks found in this matrix, this routine
%   does a more exact peak finding using parabolic fit in all direction.
%   This function is an extention of parpar_2D to 3D.
%
%     input
%         ccfmat{iPeak} : cross-correlation matrix (for all frames)
%         peaks  : peaks for all frames
%         inputs : addition, not used input
%
%     outputs
%         ccfMax  : maximumv alue of ccf
%         peaks   : peaks found after parpar
%         outputs : additonal struct as output

%  Modfications
%     12-apr-2016  GH     update ccfmat to ccfmat.data
%     11-nov-2015  GH     initial version (parpar_2D, JM)
%     11-nov-2015  GH     extention to 3D

%     11-nov-2015  JM     bug found ccfmat{iPeak} is 3D (third dim = index)
%                         inputs not used, so removed
%     28-jul-2015  JM     changed after input/output modifications
%     13-aug-2015  JM     changed, because peaks chanbed
%     07-okt-2015  JM     now row-vectors


  narginchk(2,3);
  nargoutchk(2,3);
  
  nrPeaks = length(ccfmat.data);
  ccfMax = NaN * ones(1,nrPeaks);
   
  for iPeak = 1:nrPeaks    
    %get ax,lat,ele indices of peak
    xtop = peaks.ax(iPeak);
    ytop = peaks.lat(iPeak);
    ztop = peaks.ele(iPeak);
    
    %get current max peak (non parparpar)
    ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop,ztop);
    if ccfMax(iPeak) < 0.9999

      if xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && ztop > 1 && ztop < size(ccfmat.data{iPeak},3)
        %check if it is possible to perform parabolic fitting: at least 3
        %points necessary.
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
        %x+y, z not
        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
       
      elseif xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && size(ccfmat.data{iPeak},2) < 3 && ztop > 1 && ztop < size(ccfmat.data{iPeak},3) && size(ccfmat.data{iPeak},1) >= 3 && size(ccfmat.data{iPeak},3) >= 3
          %x+z, y not

        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);     
                    
      elseif size(ccfmat.data{iPeak},1) < 3 && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && ztop > 1 && ztop < size(ccfmat.data{iPeak},3) && size(ccfmat.data{iPeak},2) >= 3 && size(ccfmat.data{iPeak},3) >= 3
        %y+z, x not
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                    
      elseif size(ccfmat.data{iPeak},1) >= 3 && xtop > 1 && xtop < size(ccfmat.data{iPeak},1) && size(ccfmat.data{iPeak},2) < 3 && size(ccfmat.data{iPeak},3) < 3
        %x, y,z not
        peaks.ax(iPeak) = (1/2*ccfmat.data{iPeak}(xtop-1,ytop,ztop)-1/2*ccfmat.data{iPeak}(xtop+1,ytop,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop)) + xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop-1,ytop,ztop)+ccfmat.data{iPeak}(xtop+1,ytop,ztop))*(ccfmat.data{iPeak}(xtop-1,ytop,ztop) - ccfmat.data{iPeak}(xtop+1,ytop,ztop))^2 +...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                    
      elseif size(ccfmat.data{iPeak},1) < 3 && size(ccfmat.data{iPeak},2) >= 3 && ytop > 1 && ytop < size(ccfmat.data{iPeak},2) && size(ccfmat.data{iPeak},3) < 3
        %y, x,z not
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop-1,ztop)-1/2*ccfmat.data{iPeak}(xtop,ytop+1,ztop))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop)) + ytop;
        peaks.ele(iPeak) = ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop-1,ztop)+ccfmat.data{iPeak}(xtop,ytop+1,ztop))*(ccfmat.data{iPeak}(xtop,ytop-1,ztop) - ccfmat.data{iPeak}(xtop,ytop+1,ztop))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);
                    
                    
       elseif size(ccfmat.data{iPeak},1) < 3 && size(ccfmat.data{iPeak},2) < 3 && size(ccfmat.data{iPeak},3) >= 3 && ztop > 1 && ztop < size(ccfmat.data{iPeak},3)
        %z, x,y not
        peaks.ax(iPeak) =  xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = (1/2*ccfmat.data{iPeak}(xtop,ytop,ztop-1)-1/2*ccfmat.data{iPeak}(xtop,ytop,ztop+1))/ ...
            (-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1)) + ztop;
        
        ccfMax(iPeak) = -1/8/(-2*ccfmat.data{iPeak}(xtop,ytop,ztop)+ccfmat.data{iPeak}(xtop,ytop,ztop-1)+ccfmat.data{iPeak}(xtop,ytop,ztop+1))*(ccfmat.data{iPeak}(xtop,ytop,ztop-1) - ccfmat.data{iPeak}(xtop,ytop,ztop+1))^2 + ...
                        ccfmat.data{iPeak}(xtop,ytop,ztop);


      else
        %z,y,x not
        peaks.ax(iPeak)  = xtop;
        peaks.lat(iPeak) = ytop;
        peaks.ele(iPeak) = ztop;
        ccfMax(iPeak) = ccfmat.data{iPeak}(xtop,ytop,ztop);
        
      end
    end    
  end
   
end
  