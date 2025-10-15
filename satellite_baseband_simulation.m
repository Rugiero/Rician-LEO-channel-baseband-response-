function [Gss taus LoSloc, signals] = satellite_baseband_simulation()
  ## The algorithm simulates the baseband equivalent signal in the LEO
  ## channel when the satellite passes over a transmitter. The LEO BS is
  ## at zenith w.r.t. its sub-satellite point. The multipath propagation is 
  ## modeled with Poisson distributed scatterers around the transmitter.
  ## At each scatterer, the signal is reflected at a random phase towards 
  ## the satellite, and the propagation delay is calculated according to
  ## the scatterer location. (The propagation between the transmitter 
  ## and each scatterer is not modeled.)
  ## 
  ## Outputs a plot of the signal and the following observables:
  ## Gss: the maximal gain of the propagated signal components at each simu-
  ##   lated time instance given in Gs(1, :), the minimum gains in Gs(2, :)
  ## taus: the maximal propagation delays at each simulated time instance
  ##       in tau(1, :), and the minimum delays in tau(2, :)
  ## LoSloc: the location of the Earth transmitter (LoS ray) at each
  ##              simulated time instance in polar coordinates
  ## signals: a vector of the received faded signal in discrete times
  
  R = 6378 * 1000; # Radius of Earth in m
  h = 200 * 1000; # Altitude of the satellite
  K = 1; # Rician parameter

  scatterers = 1000; # Number of scatterer-obstacles
  refs = refpoints(scatterers); # Initialize scatterer locations into 'refs'
  A = h; # Amplitude of the transmitted signal
  fm = 10; # Modulation frequency
  fc = 1 * 10 ^ 5; # Carrier frequency

  ## bbsignals(.) is utilized by the nested function RXbaseband(.)
  bbsignals = @(t) [A * sqrt(K / (K + 1)) .* exp(-i * 2 * pi .* fm * t)];
  ## Random phases for initial baseband transmissions
  randphases=rand(1, length(refs(1, :)) - 1);
  ## Combine the LoS and scattered components
  bbsignals = @(t) [bbsignals(t(1)) A / (sqrt(scatterers).*sqrt(1 + K))...
		    .* exp(-i * 2 * pi .* t(2: length(t))...
			    .* fm + randphases .* 2 * pi * i)];
  ## Distance to the satellites in m, utilized in RXbaseband(.)
  d = @(gamma) sqrt((cos(gamma) .* (R + h) - R) .^ 2 + ...
		    (sin(gamma) .* (R + h)) .^ 2);
  ## The angle in the antenna pattern, utilized in RXbaseband(.)
  varphi = @(gamma) acos((d(gamma) .^2 + (R + h) .^ 2 - R ^ 2)...
			 ./ (2 .* d(gamma) .* (R + h)));
  ## Gaussian antenna gain, utilized in RXbaseband(.)
  G = @(gamma) 2 .^ (-varphi(gamma) .^ 2 / deg2rad(1.6) .^ 2);
  
  N = 6000; # Number of time samples
  t = 0; # Initial time
  refs = rotateearth(-3, refs); # Initialize the Earth position
  [signal Gs] = RXbaseband(t); # Initialize the first baseband signals
  thop = 1 / 1000; # The time hop per each rotation

  taus = zeros(2, N); # Reserve the memory for the output vectors
  Gss = zeros(2, N);
  LoSloc = zeros(2, N);

  for iii = 1 : N
    ## Observe the progress
    if(mod(iii, 25000) == 0)
      iii
    end
    
    Gss(1, iii) = max(Gs); # Save the maximal gain
    Gss(2, iii) = min(Gs); # Save the minimum gain
    taus(1, iii) = max(tau); # Save the maximal propagation delay
    ## taus(2, iii) = min(tau); # Save the minimum propagation delay
    taus(2, iii) = tau(1); # Save the propagation delay of the LoS ray
    LoSloc(:, iii) = refs(:, 1); # Save the transmitter location
    
    refs = rotateearth(thop, refs); # Rotate Earth.
    [signal Gs] = RXbaseband(t); # New received baseband signal
    signals(iii) = signal;
    t = t + thop; # Next time instance
  end

  ## Plot the resulting signal
  figure;
  hold on;
  plot(linspace(0, t, N), real(signals), 'linewidth', 1);
  axis([[0, t], [-1, 1]]);
  xlabel('time (s)', 'fontname', 'DejaVu Serif', 'fontsize', 12);
  legend('Sinusoidal baseband signal', 'fontname', 'DejaVu Serif',...
	 'interpreter', 'tex', 'fontsize', 12);
  grid on;
  ylabel('Amplitude (non-dimensional)', 'fontname',...
	 'DejaVu Serif', 'interpreter', 'tex', 'fontsize', 12);
  set(gca,'fontsize',12);
  title('LoS channel, f_m = 10 Hz, f_c = 0.1 GHz',...
	'fontname', 'DejaVu Serif', 'interpreter', 'tex', 'fontsize', 12);
  hold off;

  ## In the following, define two nested functions for the derivation of the   ## baseband signal at each time instance, and other for the Earth rotation

  ## The received baseband signal at time t
  function [signal Gs] =RXbaseband(t)
    c = 299792458; # Speed of light
    if(!isempty(refs))
      a = 1 ./ d(refs(1, :)); # Spatial amplitude path losses
      Gs= G(refs(1, :)); # Gaussian antenna pattern
      ## Gs=ones(1, length(refs(1, :))); # Omnidirectional antenna pattern
      tau  = d(refs(1, :)) / c; # Propagation delays
    else
      a = [0];
      tau = [0];
    end
    ## Received individual signals
    ab = Gs .* a .* exp(-i * 2 * pi .* tau .* fc);
    signal = sum(ab .* bbsignals(t - tau)); #Aggregate received signal
  endfunction

  ## Rotates the Earth (moves the satellite) and the scatterer locations in    ## 'refs' over the period thop, which can be negative
  
  function newrefs = rotateearth(thop, refs)   
    GM = 3.986 * 10 ^ 14; # Gravitational constant
    orbitalspeed = sqrt(GM / (h + R)); # Satellite speed in m/s
    angularspeed = orbitalspeed / (h + R); # Angular speed of the satellite
    rotation = angularspeed * thop; # Rotation of Earth
    ## Transform the polar coordinates to euclidean coordinates
    eucpos  = pol2euc(refs);
    ## Rotation about x-axis
    newpos = [[1 0 0]; [0 cos(rotation) -sin(rotation)];...
	      [0 sin(rotation) cos(rotation)]] * eucpos;
    newrefs = euc2pol(newpos); # Back to the polarcoordinates
  endfunction

endfunction


## Returns a table of numbPoints Poisson points around the north pole
function refs = refpoints(numbPoints)
  refs = [0; 0]; # LOS component location.
  ## Generate the Poisson distributed random obstacle locations
  ## Change yMin to control the area width of the scatterers
  ##  yMin = 1 - 0.00035; yMax = 1; # -10 dB footprint for h = 2000 km
  yMin = 1 - 0.00000035; yMax = 1; # -3 dB footprint for h = 200 km
  xMin = -pi; xMax = pi;
  xDelta = xMax - xMin; yDelta = yMax - yMin; # Rectangle dimensions
  ## Pick points from uniform distribution
  x = xDelta * (rand(numbPoints, 1)) + xMin;    
  y = yDelta * (rand(numbPoints, 1)) + yMin;    
  ## Map referencepoints to spherical coordinates
  refs = [refs [pi / 2 - asin(y)'; x']]; 
end

function p  = euc2pol(e) # Euclidean 3D coordinates to polar
  R = 6378 * 1000; # Radius of Earth
  p = [acos(e(3, :) ./ R); (e(2, :) >= 0) .* atan2(e(2, :), e(1, :)) ...
			   + (e(2, :) < 0) .* (atan2(e(2, :), e(1, :)) ...
					       + 2 * pi)]; 
end

function e = pol2euc(p) # Polar coordinates to 3D Euclidean
  R= 6378 * 1000; # Radius of Earth.
  e = [R * cos(p(2, :)) .* sin(p(1, :));...
       R * sin(p(2, :)) .* sin(p(1, :));...
       R * cos(p(1, :))];
end
