pro irxbeta, ir_struc, irflux_sub, irfluxerr_sub, $
             opt_struc, optflux_sub, optfluxerr_sub, optwave, $
             z_struc, z_sub, ztype_sub, list, $
             index, index_err, f1600, f1600_err, l1600, l1600_err, lir, $
             lir_err, method, pred_flux_870 = pred_flux_870

; INPUTS

; ir_struc = the IR data structure, containing the following fluxes in
; uJy: 1) MIPS 24, 2) MIPS 70, 3) PACS 70, 4) PACS 100, 5) PACS 160,
; 6) SPIRE 250, 7) SPIRE 350, 8) SPIRE 500, 9) ALMA 870

; irflux_sub = the subscripts from TAG_NAMES of the fluxes for IR SED
; fitting

; irfluxerr_sub = for flux errors

; opt_struc = the optical data structure, containing RA, Dec, at least
; 2 fluxes in uJy

; optflux_sub = the subscripts from TAG_NAMES of the fluxes for beta
; fitting, etc.

; optfluxerr_sub = for flux errors

; optwave = the wavelengths in micron corresponding to the optical fluxes

; z_struc = the redshift data structure, ELEMENTS MUST BE SAME AS IN
; OPTICAL STRUCTURE

; z_sub = the subscripts from TAG_NAMES of the redshifts

; ztype_sub = the subscripts from TAG_NAMES of the redshift types

; list = the subscripts in the optical structure matched to the ir structure

; OUTPUTS

; index

; index error

; f1600

; f1600 error

; l1600 in Lsun

; l1600 error in Lsun

; LIR in Lsun

; method = 1 if used >3 bands in fit for index; 2 if used 2 bands in
; fit for index; 3 if used 1 band + limit in fit for index (resulting
; in a lower limit for index); 4 if used 1 band + limit in fit for
; index (resulting in an upper limit for index)

forward_function ce01_general

!path = lir_path+!path

index = fltarr(n_elements(list))+!values.f_nan
index_err = fltarr(n_elements(list))+!values.f_nan
normalization = dblarr(n_elements(list))+!values.f_nan
;normalization_err = fltarr(n_elements(list))+!values.f_nan
method = fltarr(n_elements(list))
f1600 = fltarr(n_elements(list))+!values.f_nan
f1600_err = fltarr(n_elements(list))+!values.f_nan
l1600 = fltarr(n_elements(list))+!values.f_nan
l1600_err = fltarr(n_elements(list))+!values.f_nan
lir = fltarr(n_elements(list))+!values.f_nan
lir_err = fltarr(n_elements(list))+!values.f_nan
pred_flux_870 = fltarr(n_elements(list))+!values.f_nan

for i=0,n_elements(list)-1 do begin

    if finite(list[i]) eq 1 then begin

        allz = dblarr(n_elements(z_sub))
        allztype = strarr(n_elements(z_sub))

        for m=0,n_elements(z_sub)-1 do begin
            
            allz[m] = z_struc.(z_sub[m])[list[i]]
            allztype[m] = z_struc.(ztype_sub[m])[list[i]]
            
        endfor
        
        zspec = where(strmatch(allztype, 'S') eq 1, numzspec)
        zphot = where(strmatch(allztype, 'P') eq 1, numzphot)
        
        if numzspec eq 1 then redshift = allz[zspec[0]]
        if numzspec eq 0 and numzphot eq 1 then redshift = allz[zphot[0]]

        if n_elements(redshift) eq 1 then begin

           if redshift ne 0 then begin

              optflux = dblarr(n_elements(optflux_sub))
              optfluxerr = dblarr(n_elements(optflux_sub))
              wave = double(optwave*10000.)

              for j=0,n_elements(optflux_sub)-1 do begin

                 optflux[j] = opt_struc.(optflux_sub[j])[list[i]]
                 optfluxerr[j] = opt_struc.(optfluxerr_sub[j])[list[i]]

              endfor

              order = sort(wave)
              wave = wave[order]
              optflux = optflux[order]
              optfluxerr = optfluxerr[order]

              range = where(wave/(1.+redshift) lt 2630. and $
                            wave/(1.+redshift) gt 1250. and $
                            optfluxerr gt 0, numrange, $
                            complement = out_range, ncomplement = numout_range)

              positive = where(optflux[range] gt 0, numrange_pos)

;             estimate with 2 good data points

              if numrange ge 2 and numrange_pos ge 2 then begin

                 if numout_range gt 0 then $
                    remove, out_range, optflux, optfluxerr, wave

                 subrange_pos = where(optflux gt 0.)
                 min_subrange_pos = min(subrange_pos)
                 max_subrange_pos = max(subrange_pos)

                 wave_norm = wave[min_subrange_pos]

;                remove scale from fitting

                 wave = wave/wave_norm

                 optflux_norm = optflux[min_subrange_pos]

                 optflux_lambda = optflux/(wave^2.)/optflux_norm

                 optfluxerr_lambda = optfluxerr/(wave^2.)/optflux_norm

;                initial guess

                 p0 = [1., alog10(1./optflux_lambda[max_subrange_pos])/$
                           alog10(1./wave[max_subrange_pos])]

;                parameter constraints

                 pinfo = replicate({limited: [1,0], limits: [0.d,0]}, 2)
                 pinfo[1].limited[0] = 0

                 args = {wave: wave, flux_lambda: optflux_lambda, $
                         err_lambda: optfluxerr_lambda}

;                non-linear regression

                 fit = mpfit('beta', p0, functargs = args, $
                             bestnorm = bestnorm, dof = dof, $
                             niter = niter, maxiter = 200, $
                             parinfo = pinfo, covar = covar, /quiet)

                 c = 3d10 ; in cm*Hz

                 common_factor = c*1d-21*optflux_norm/(wave_norm^(fit[1]+2.))

                 normalization[i] = fit[0]*common_factor

;                 the following formula is incorrect but included here
;                 for historical purposes
;                 normalization_err[i] = sqrt((common_factor*param_err[0])^2.+$
;                   (common_factor*fit[0]*alog(wave_norm)*param_err[1])^2.+$
;                   2.*fit[0]*common_factor^2.*alog(wave_norm)*covar[0,1])

                 index[i] = fit[1]

                 if numrange gt 2 then begin
                 
                   covar = covar*bestnorm/dof
                   add_ssr = total((normalization[i]*(wave*wave_norm)^index[i]-$
                            optflux_lambda*common_factor*wave_norm^index[i])^2.)
                   
                 endif else if numrange eq 2 then begin
                 
                   add_ssr = 0d
                   
                 endif

                 index_err[i] = sqrt(covar[1,1])

                 f1600[i] = normalization[i]*(1600.*(1.+redshift))^(index[i])

                 common_factor_1600 = c*1d-21*optflux_norm*$
                   (1600.*(1.+redshift))^(fit[1])/(wave_norm^(fit[1]+2.))

                 f1600_err[i] = sqrt(covar[0,0]*(common_factor_1600)^2.+$
                   covar[1,1]*(fit[0]*common_factor_1600*$
                   alog(1600.*(1.+redshift)/wave_norm))^2.+$
                   2.*fit[0]*common_factor_1600^2.*$
                   alog(1600.*(1.+redshift)/wave_norm)*$
                   covar[0,1]+add_ssr)

;                   sqrt(1.+1./n_elements(wave)+$
;                   (1600.*(1.+redshift)/wave_norm-mean(wave))^2/$
;                   total((wave-mean(wave))^2))

;                  prediction interval term for linear regression

                 dist = dlum_all(redshift)

                 l1600[i] = (1600.*(1.+redshift))*f1600[i]*4*!pi*$
                            (dist*3.0856d24)^2/3.826d33

                 l1600_err[i] = (1600.*(1.+redshift))*f1600_err[i]*4*!pi*$
                                (dist*3.0856d24)^2/3.826d33

;                for inspection of cases when fitting doesn't converge

                 if total(covar eq 0) ne 0 then stop
                 if niter eq 200 then stop

                 if numrange gt 2 then method[i] = 1 else method[i] = 2

              endif

;             estimate with limits, lower limit for index

              if numrange ge 2 and numrange_pos lt 2 then begin

                 subrange1 = where(wave/(1.+redshift) lt 2630. and $
                                   wave/(1.+redshift) gt 1250. and $
                                   optflux/optfluxerr ge 3., numsubrange1)

                 subrange2 = where(wave/(1.+redshift) lt 2630. and $
                                   wave/(1.+redshift) gt 1250. and $
                                   optfluxerr gt 0.)

                 if numsubrange1 gt 0 and min(subrange2) ne subrange1 $
                 and max(subrange2) eq subrange1 then begin

                    min_subrange = min(subrange2)
                    max_subrange = subrange1

                    optflux_new = [3.*optfluxerr[min_subrange], $
                                   optflux[max_subrange]]

                    optflux = optflux_new

                    wave_new = [wave[min_subrange], wave[max_subrange]]

                    wave = wave_new

                    c = 3d18

                    optflux_lambda = c/(wave^2)*optflux*1.d-29

;                   solution is analytic in this case

                    index[i] = alog10(optflux_lambda[0]/optflux_lambda[1])/alog10(wave[0]/wave[1])
                    
                    f1600[i] = optflux_lambda[0]*(1600.*(1.+redshift)/wave[0])^(index[i])

                    dist = dlum_all(redshift)

                    l1600[i] = (1600.*(1.+redshift))*f1600[i]*4*!pi*(dist*3.0856d24)^2/3.826d33

                    method[i] = 3
           
                 endif

;                estimate with limits, upper limit for index

                 if numsubrange1 gt 0 and max(subrange2) ne subrange1 $
                 and min(subrange2) eq subrange1 then begin

                    min_subrange = subrange1
                    max_subrange = max(subrange2)

                    optflux_new = [optflux[min_subrange], $
                                   3.*optfluxerr[max_subrange]]

                    optflux = optflux_new

                    wave_new = [wave[min_subrange], wave[max_subrange]]

                    wave = wave_new

                    c = 3d18

                    optflux_lambda = c/(wave^2)*optflux*1.d-29

                    index[i] = alog10(optflux_lambda[0]/optflux_lambda[1])/alog10(wave[0]/wave[1])
                    
                    f1600[i] = optflux_lambda[0]*(1600.*(1.+redshift)/wave[0])^(index[i])

                    dist = dlum_all(redshift)

                    l1600[i] = (1600.*(1.+redshift))*f1600[i]*4*!pi*(dist*3.0856d24)^2/3.826d33

                    method[i] = 4
           
                 endif

              endif

              irflux = fltarr(n_elements(irflux_sub))
              irfluxerr = fltarr(n_elements(irflux_sub))

              for k=0,n_elements(irflux_sub)-1 do begin

                 irflux[k] = ir_struc.(irflux_sub[k])[i]
                 irfluxerr[k] = ir_struc.(irfluxerr_sub[k])[i]

              endfor

;             get lir

              lir[i] = ce01_general(irflux, irfluxerr, redshift, $
                                    lir_err = lir_err1, $
                                    pred_fluxes = pred_fluxes_1)
              lir_err[i] = lir_err1

              pred_flux_870[i] = pred_fluxes_1[n_elements(pred_fluxes_1)-1]

           endif

        endif

;       make redshift undefined again

        redshiftgoaway = size(temporary(redshift))

     endif


 endfor

end
