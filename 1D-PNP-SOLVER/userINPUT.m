phi0     = 5    ;% applied potential on the electrode at z = 0 (kBT/e)
ac       = 1    ;% [0 for dc, 1 for ac]
freq     = 100e3;% frequency (Hz) [will be ignored for dc]
x_M      = 1e-3 ;% electrolyte strength (M)
D1       = 1e-9 ;% cation diffusivity (m^2/s)
D2       = 2e-9 ;% anion diffusivity (m^2/s)
q1       = +1   ;% + ion charge number (-)
q2       = -1   ;% - ion charge number (-)
epsinf   = 78.36;% dielectric constant (-)
H        = 100e-9;% gap (m)
h_S      = 0    ;% thickness of the Stern layer (m)
epsRatio = 0.1  ;% Stern to diffuse layer permittivity ratio (-)