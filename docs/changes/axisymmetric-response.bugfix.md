Dual-reflector optical-model compilation now reads primary reflectivity from the mirror record rather than its geometric surface prescription. Complete wavelength/incidence-angle response grids are preserved and validated for native transport; legacy CSV model export rejects angle-dependent response instead of discarding it.

Transform the secondary polynomial from sim_telarray's reversed local mirror
frame into telescope coordinates. Nominal SST and SCT smoke rays now reach the
focal surface. Preserve measured reflectivity metadata alongside the response;
these runs do not constitute production qualification.
