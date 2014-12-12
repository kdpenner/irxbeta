function elb11_general, fluxes, flux_errs, redshift, lir_err = elb11_lir_err

; fluxes must have entries in the following order:
; mips 24, mips70, pacs 70, pacs 100, pacs 160, spire 250, spire 350,
; spire 500.  Use NaN if the flux and its error are not available.

readcol, 'mips_24_response.dat', wave_mips24, resp_mips24, $
format = '(d,d)', /silent
readcol, 'mips_70_response.dat', wave_mips70, resp_mips70, $
format = '(d,d)', /silent
readcol, 'pacs_70_response.dat', wave_pacs70, resp_pacs70, $
format = '(d,d)', /silent
readcol, 'pacs_100_response.dat', wave_pacs100, resp_pacs100, $
format = '(d,d)', /silent
readcol, 'pacs_160_response.dat', wave_pacs160, resp_pacs160, $
format = '(d,d)', /silent
readcol, 'spire_250_response.dat', wave_spire250, resp_spire250, $
format = '(d,d)', /silent
readcol, 'spire_350_response.dat', wave_spire350, resp_spire350, $
format = '(d,d)', /silent
readcol, 'spire_500_response.dat', wave_spire500, resp_spire500, $
format = '(d,d)', /silent

H_0 = 70.
o_m = 0.3
o_l = 0.7
c_light = 2.99793d14

distance = dlum_all(redshift, H0=H0, o_m=o_m, o_l=o_l)

readcol, 'sed_elbaz2011_ms.nuLnu', lambda_ms, nuLnuinLsun_ms, $
format = '(d,d)', /silent
readcol, 'sed_elbaz2011_sb.nuLnu', lambda_sb, nuLnuinLsun_sb, $
format = '(d,d)', /silent

lambda = lambda_ms

nuLnuinLsun = transpose([[nuLnuinLsun_ms], [nuLnuinLsun_sb]])

lir = [1.d11, 1.d11]

lambda_obs = lambda*(1.+redshift)

nuSnu = nuLnuinLsun/(4*!dpi*(distance*3.0856d22)^2/3.826d26)

num_models = (size(nuLnuinLsun, /dimen))[0]

Snu_mips24 = dblarr(num_models)
Snu_mips70 = dblarr(num_models)
Snu_pacs70 = dblarr(num_models)
Snu_pacs100 = dblarr(num_models)
Snu_pacs160 = dblarr(num_models)
Snu_spire250 = dblarr(num_models)
Snu_spire350 = dblarr(num_models)
Snu_spire500 = dblarr(num_models)

;FUCK YOU INT_TABULATED

bbody_ref_mips24 = planck(wave_mips24*1.d4, 10000.)
bbody_ref_mips70 = planck(wave_mips70*1.d4, 10000.)

norm_mips24 = tsum(wave_mips24, resp_mips24*bbody_ref_mips24)/$
              planck(23.675*1.d4, 10000.)
norm_mips70 = tsum(wave_mips70, resp_mips70*bbody_ref_mips70)/$
              planck(71.440*1.d4, 10000.)
norm_pacs70 = 1.d-32*tsum(c_light/wave_pacs70, resp_pacs70*c_light/wave_pacs70)
norm_pacs100 = 1.d-32*tsum(c_light/wave_pacs100, $
                           resp_pacs100*c_light/wave_pacs100)
norm_pacs160 = 1.d-32*tsum(c_light/wave_pacs160, $
                           resp_pacs160*c_light/wave_pacs160)
norm_spire250 = 1.d-32*tsum(c_light/wave_spire250, $
                            resp_spire250*c_light/wave_spire250)
norm_spire350 = 1.d-32*tsum(c_light/wave_spire350, $
                            resp_spire350*c_light/wave_spire350)
norm_spire500 = 1.d-32*tsum(c_light/wave_spire500, $
                            resp_spire500*c_light/wave_spire500)

for i=0,num_models-1 do begin

   nuSnu_interp_mips24 = interpol(nuSnu[i,*], lambda_obs, wave_mips24)
   Slambda_interp_mips24 = nuSnu_interp_mips24/(10.*wave_mips24)
   Slambda_mips24 = tsum(wave_mips24, resp_mips24*Slambda_interp_mips24)/$
                    norm_mips24
   Snu_mips24[i] = 1.d29*1.d4/(c_light/23.675^2)*Slambda_mips24


   nuSnu_interp_mips70 = interpol(nuSnu[i,*], lambda_obs, wave_mips70)
   Slambda_interp_mips70 = nuSnu_interp_mips70/(10.*wave_mips70)
   Slambda_mips70 = tsum(wave_mips70, resp_mips70*Slambda_interp_mips70)/$
                    norm_mips70
   Snu_mips70[i] = 1.d29*1.d4/(c_light/71.440^2)*Slambda_mips70


   nuSnu_interp_pacs70 = interpol(nuSnu[i,*], lambda_obs, wave_pacs70)
   Snu_pacs70[i] = tsum(c_light/wave_pacs70, resp_pacs70*nuSnu_interp_pacs70)/$
                   norm_pacs70


   nuSnu_interp_pacs100 = interpol(nuSnu[i,*], lambda_obs, wave_pacs100)
   Snu_pacs100[i] = tsum(c_light/wave_pacs100, $
                         resp_pacs100*nuSnu_interp_pacs100)/norm_pacs100


   nuSnu_interp_pacs160 = interpol(nuSnu[i,*], lambda_obs, wave_pacs160)
   Snu_pacs160[i] = tsum(c_light/wave_pacs160, $
                         resp_pacs160*nuSnu_interp_pacs160)/norm_pacs160


   nuSnu_interp_spire250 = interpol(nuSnu[i,*], lambda_obs, wave_spire250)
   Snu_spire250[i] = tsum(c_light/wave_spire250, $
                          resp_spire250*nuSnu_interp_spire250)/norm_spire250


   nuSnu_interp_spire350 = interpol(nuSnu[i,*], lambda_obs, wave_spire350)
   Snu_spire350[i] = tsum(c_light/wave_spire350, $
                          resp_spire350*nuSnu_interp_spire350)/norm_spire350


   nuSnu_interp_spire500 = interpol(nuSnu[i,*], lambda_obs, wave_spire500)
   Snu_spire500[i] = tsum(c_light/wave_spire500, $
                          resp_spire500*nuSnu_interp_spire500)/norm_spire500

endfor

model_fluxes = [[Snu_mips24], [Snu_mips70], [Snu_pacs70], [Snu_pacs100], $
                [Snu_pacs160], [Snu_spire250], [Snu_spire350], [Snu_spire500]]

good = where(finite(fluxes) eq 1, good_count)

if good_count gt 1 then begin

   chisq_models = dblarr(num_models)+!values.d_nan
   factor_models = dblarr(num_models)
   factor_err_models = dblarr(num_models)

   for j=0,num_models-1 do begin

      pstart = [median(fluxes/model_fluxes[j,*])] > [0.01d]

      pinfo = {limited: [1,0], limits: [0.d,0]}

      factor = mpfitfun('template_renorm', model_fluxes[j,*], fluxes, $
                        flux_errs, pstart, parinfo = pinfo, bestnorm = chisq, $
                        perror = factor_err, /nan, /quiet)

      chisq_models[j] = chisq
      factor_models[j] = factor
      factor_err_models[j] = factor_err*sqrt(chisq/(good_count-1))

   endfor

   best_model = where(chisq_models eq min(chisq_models))
   elb11_lir = factor_models[best_model]*lir[best_model]
   elb11_lir_err = factor_err_models[best_model]*lir[best_model]

   return, elb11_lir

endif

if good_count eq 1 then begin

   elb11_lir = fluxes[good[0]]/model_fluxes[0,good]*lir[0]
   elb11_lir_err = flux_errs[good[0]]/model_fluxes[good[0]]*lir[0]

   return, elb11_lir

endif

if good_count eq 0 then begin

   elb11_lir_err = 0.
   return, 0.

endif


end
