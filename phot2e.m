

function E = phot2e(lambda,molphotons)
%E = phot2e(lambda,molphotons) calculates the E Joules of energy corresponding to the number of moles of photons of wavelength lambda (m)
A         = 6.02214E23; % [mol-1]       Constant of Avogadro
h         = 6.6262E-34; % [J s]         Planck's constant
c         = 299792458;  % [m s-1]       Speed of light
%calculates the energy content (J) of 1 photon of wavelength lambda (m)
e           =h*c./lambda;   % [J]           energy of 1 photon
photons  = A.*molphotons;
E     = photons.*e*1000;

end
