function [rada, radp] = genRadiationPattern4(beamwidth)

  Nt=100;
  deltat=1;
  Lt=0.2;
  step = 0.001;
  omegat=-pi:step:pi + step;
  gainfunction = ones(size(omegat));

  mingain = 0.0001;

  for ii=1:length(omegat)
	   %if (omegat(ii) == 0)
	    %	gainfunction(ii) = 1;
	     %else
		     %gainfunction(ii)= (1/Nt) * exp(1i*pi*deltat* omegat(ii)* (Nt-1)) * (sin(pi*Lt*omegat(ii))/sin(pi*Lt*omegat(ii)*Nt^-1).^1);
		     %if (abs(omegat(ii))  > beamwidth)
	       %   gainfunction(ii) = 0;
	       %else
	      %%%%    gainfunction(ii) = sqrt(sin(omegat(ii) - pi/2)); %2 - (2/beamwidth)*(abs(omegat(ii)) );
        gainfunction(ii) = (sin(omegat(ii) - pi/2)).^beamwidth;
			      %gainfunction(ii) = -omegat(ii)^2 + (beamwidth/2)^2;
	      %end	%end


  end
  % have a minimum gain level
  gainfunction = gainfunction + 0.001;
  % gain cannot be negative
  gainfunction(find(gainfunction < mingain)) = mingain;

  gainfunction = gainfunction / max(gainfunction);
  gainfunction(find(omegat == 0)) = 1;

  rada = omegat;
  radp = gainfunction;
end
