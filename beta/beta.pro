function beta, p, wave = wave, flux_lambda = flux_lambda, $
               err_lambda = err_lambda

; wave is in Angstrom

model = p[0]*(wave)^(p[1])

err = err_lambda

obs = flux_lambda

return, (obs - model)/err


end
