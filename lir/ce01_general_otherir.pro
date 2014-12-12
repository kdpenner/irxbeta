function ce01_general_otherir, fluxes, flux_errs, redshift, $
                               lir_err = ce_lir_err

; fluxes must have entries in the following order:
; 1) AKARI 9, 2) WISE 12, 3) WISE 22, 4) IRAS 25, 5) IRAS 60,
; 6) AKARI 65, 7) AKARI 90, 8) IRAS 100, 9) AKARI 140
; Use NaN if the flux and its error are not available.

readcol, 'irc_9_response.dat', wave_irc9, resp_irc9, $
format = '(d,d)', /silent
readcol, 'wise_12_response.dat', wave_wise12, resp_wise12, $
format = '(d,d)', /silent
readcol, 'wise_22_response.dat', wave_wise22, resp_wise22, $
format = '(d,d)', /silent
readcol, 'iras_25_response.dat', wave_iras25, $
resp_iras25, format = '(d,d)', /silent
readcol, 'iras_60_response.dat', wave_iras60, $
resp_iras60, format = '(d,d)', /silent
readcol, 'fis_65_response.dat', wave_fis65, $
resp_fis65, format = '(d,d)', /silent
readcol, 'fis_90_response.dat', wave_fis90, $
resp_fis90, format = '(d,d)', /silent
readcol, 'iras_100_response.dat', wave_iras100, $
resp_iras100, format = '(d,d)', /silent
readcol, 'fis_140_response.dat', wave_fis140, $
resp_fis140, format = '(d,d)', /silent

H_0 = 70.
o_m = 0.3
o_l = 0.7
c_light = 2.99793d14

distance = dlum_all(redshift, H0=H0, o_m=o_m, o_l=o_l)

restore, 'chary_elbaz.save'

lambda_obs = lambda*(1.+redshift)

nuSnu = nuLnuinLsun/(4*!dpi*(distance*3.0856d22)^2/3.826d26)

num_models = (size(nuLnuinLsun, /dimen))[0]

Snu_irc9 = dblarr(num_models)
Snu_wise12 = dblarr(num_models)
Snu_wise22 = dblarr(num_models)
Snu_iras25 = dblarr(num_models)
Snu_iras60 = dblarr(num_models)
Snu_fis65 = dblarr(num_models)
Snu_fis90 = dblarr(num_models)
Snu_iras100 = dblarr(num_models)
Snu_fis140 = dblarr(num_models)

;FUCK YOU INT_TABULATED

norm_irc9 =  1.d-32*tsum(c_light/wave_irc9, $
                          resp_irc9*c_light/wave_irc9)
norm_wise12 = 1.d-32*tsum(c_light/wave_wise12, $
                          resp_wise12*c_light/wave_wise12)
norm_wise22 = 1.d-32*tsum(c_light/wave_wise22, $
                          resp_wise22*c_light/wave_wise22)
norm_iras25 = 1.d-32*tsum(c_light/wave_iras25, $
                           resp_iras25*c_light/wave_iras25)
norm_iras60 = 1.d-32*tsum(c_light/wave_iras60, $
                           resp_iras60*c_light/wave_iras60)
norm_fis65 = 1.d-32*tsum(c_light/wave_fis65, $
                            resp_fis65*c_light/wave_fis65)
norm_fis90 = 1.d-32*tsum(c_light/wave_fis90, $
                            resp_fis90*c_light/wave_fis90)
norm_iras100 = 1.d-32*tsum(c_light/wave_iras100, $
                            resp_iras100*c_light/wave_iras100)
norm_fis140 = 1.d-32*tsum(c_light/wave_fis140, $
                           resp_fis140*c_light/wave_fis140)


for i=0,num_models-1 do begin

   nuSnu_interp_irc9 = interpol(nuSnu[i,*], lambda_obs, wave_irc9)
   Snu_irc9[i] = tsum(c_light/wave_irc9, $
                        resp_irc9*nuSnu_interp_irc9)/norm_irc9

   nuSnu_interp_wise12 = interpol(nuSnu[i,*], lambda_obs, wave_wise12)
   Snu_wise12[i] = tsum(c_light/wave_wise12, $
                        resp_wise12*nuSnu_interp_wise12)/norm_wise12

   nuSnu_interp_wise22 = interpol(nuSnu[i,*], lambda_obs, wave_wise22)
   Snu_wise22[i] = tsum(c_light/wave_wise22, $
                        resp_wise22*nuSnu_interp_wise22)/norm_wise22

   nuSnu_interp_iras25 = interpol(nuSnu[i,*], lambda_obs, wave_iras25)
   Snu_iras25[i] = tsum(c_light/wave_iras25, $
                         resp_iras25*nuSnu_interp_iras25)/norm_iras25

   nuSnu_interp_iras60 = interpol(nuSnu[i,*], lambda_obs, wave_iras60)
   Snu_iras60[i] = tsum(c_light/wave_iras60, $
                         resp_iras60*nuSnu_interp_iras60)/norm_iras60

   nuSnu_interp_fis65 = interpol(nuSnu[i,*], lambda_obs, wave_fis65)
   Snu_fis65[i] = tsum(c_light/wave_fis65, $
                          resp_fis65*nuSnu_interp_fis65)/norm_fis65

   nuSnu_interp_fis90 = interpol(nuSnu[i,*], lambda_obs, wave_fis90)
   Snu_fis90[i] = tsum(c_light/wave_fis90, $
                          resp_fis90*nuSnu_interp_fis90)/norm_fis90

   nuSnu_interp_iras100 = interpol(nuSnu[i,*], lambda_obs, wave_iras100)
   Snu_iras100[i] = tsum(c_light/wave_iras100, $
                          resp_iras100*nuSnu_interp_iras100)/norm_iras100

   nuSnu_interp_fis140 = interpol(nuSnu[i,*], lambda_obs, wave_fis140)
   Snu_fis140[i] = tsum(c_light/wave_fis140, $
                         resp_fis140*nuSnu_interp_fis140)/norm_fis140


endfor

model_fluxes = [[Snu_irc9], [Snu_wise12], [Snu_wise22], [Snu_iras25], $
                [Snu_iras60], [Snu_fis65], [Snu_fis90], $
                [Snu_iras100], [Snu_fis140]]

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
   ce_lir = factor_models[best_model]*lir[best_model]
   ce_lir_err = factor_err_models[best_model]*lir[best_model]

   pred_fluxes = (factor_models[best_model])[0]*model_fluxes[best_model,*]

   return, ce_lir

endif

if good_count eq 1 then begin

   minval = min(abs(fluxes[good[0]]-model_fluxes[*,good]), best_model)
   ce_lir = fluxes[good[0]]/model_fluxes[best_model,good]*lir[best_model]
   ce_lir_err = flux_errs[good[0]]/model_fluxes[best_model,good]*lir[best_model]

   pred_fluxes = (fluxes[good[0]]/model_fluxes[best_model,good])[0]*model_fluxes[best_model,*]

   return, ce_lir

endif

if good_count eq 0 then begin

   ce_lir_err = 0.

   pred_fluxes = fltarr((size(model_fluxes, /dimen))[1])

   return, 0.

endif


end
