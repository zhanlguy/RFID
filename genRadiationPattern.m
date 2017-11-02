function [rada, radp] = genRadiationPattern(beamwidth)

  Nt=100;
  deltat=1;
  Lt=0.2;
  omegat=-2:0.001:2;
  gainfunction = ones(size(omegat));
  for ii=1:length(omegat)
	   %if (omegat(ii) == 0)
	    %	gainfunction(ii) = 1;
	     %else
		     %gainfunction(ii)= (1/Nt) * exp(1i*pi*deltat* omegat(ii)* (Nt-1)) * (sin(pi*Lt*omegat(ii))/sin(pi*Lt*omegat(ii)*Nt^-1).^1);
		     if (abs(omegat(ii))  > beamwidth/2)
	          gainfunction(ii) = 0;
	       else
	          %gainfunction(ii) = 1 - (2/beamwidth)*(abs(omegat(ii)) );
			      gainfunction(ii) = -omegat(ii)^2 + (beamwidth/2)^2;
	      end	%end
  end
  gainfunction = gainfunction / max(gainfunction);
  gainfunction(find(omegat == 0)) = 1;

  rada = omegat;
  radp = gainfunction;
end
