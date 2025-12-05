program define MJleeIMEEME_BASE2, rclass
    version 14.0
    
    syntax varlist(min=3 max=3) ///
        , ///
        SELection_int(varlist) /// interaction vars in selection eq
        SELection_cov(varlist) /// covariates in selection eq
        OUTcome_int(varlist)   /// interaction vars in outcome eq
        OUTcome_cov(varlist)   /// covariates in outcome eq
        
    tokenize `varlist'
    local resvarb   `1'
    local treatvarb `2'
    local partvarb  `3'
    
    qui {
        tempvar touse d q y ad adwW wpro_2 dwpro wpro prob dens adnadwW
        tempvar cd_t cd_i cd cd2 cd3 lam2 lam3 lam4 coveu g1 g2 g3
        
        mark `touse'
        
        // Generate basic variables
        gen `d' = `treatvarb' if `touse'
        gen `q' = `partvarb' if `touse'
        gen `y' = `resvarb' if `touse'
        
        // Generate interaction terms for selection equation
        local dinteract ""
        foreach var of local selection_int {
            tempvar dvar_sel
            gen `dvar_sel' = `var' * `d' if `touse'
            local dinteract "`dinteract' `dvar_sel'"
            local dvar_map_sel_`var' "`dvar_sel'"
        }
        
        // Create lists for selection equation
        local dw "`d' `selection_cov' `dinteract'"
        
        // Generate interaction terms for outcome equation
        local doutinteract ""
        foreach var of local outcome_int {
            tempvar dvar_out
            gen `dvar_out' = `var' * `d' if `touse'
            local doutinteract "`doutinteract' `dvar_out'"
            local dvar_map_out_`var' "`dvar_out'"
        }
        
        local dx "`d' `outcome_cov' `doutinteract'"
        
        /* Probit */
        probit `q' `dw'
        gen `ad' = _b[`d']
        gen `adwW' = 0
        foreach var of local selection_int {
            replace `adwW' = `adwW' + _b[`dvar_map_sel_`var''] * `var'
        }
        tempvar wpro_xb
        predict `wpro_xb', xb
        gen `dwpro' = `wpro_xb'
        gen `wpro' = `wpro_xb' - `d'*_b[`d'] - `d'*`adwW'
        gen `prob' = normal(`wpro')
        gen `dens' = normalden(`wpro')
        gen `adnadwW' = `ad' + `adwW'

        /* Tobit model */
        tobit `y' `dx', ll(0)
        gen `cd_t' = _b[`d']
        foreach var of local outcome_int {
            replace `cd_t' = `cd_t' + _b[`dvar_map_out_`var''] * `var'
        }
        tempname tob_mean
        egen `tob_mean' = mean(`cd_t'*`prob')
        scalar tob = `tob_mean'

        /* Independence model */
        reg `y' `dx' if `q'==1
        gen `cd_i' = _b[`d']
        foreach var of local outcome_int {
            replace `cd_i' = `cd_i' + _b[`dvar_map_out_`var''] * `var'
        }
        tempvar yhat_i XB_i
        predict `yhat_i', xb
        gen `XB_i' = `yhat_i' - `cd_i'*`d'

        tempname imeind_mean emeind_mean
        egen `imeind_mean' = mean(`cd_i'*`prob')
        scalar imeind = `imeind_mean'
        egen `emeind_mean' = mean(`adnadwW'*`dens'*`XB_i')
        scalar emeind = `emeind_mean'
        scalar totind = imeind + emeind

        /* Normality model (Heckit) */
        gen `lam2' = normalden(`dwpro')/normal(`dwpro')
        reg `y' `dx' `lam2' if `q'==1
        gen `cd' = _b[`d']
        foreach var of local outcome_int {
            replace `cd' = `cd' + _b[`dvar_map_out_`var''] * `var'
        }
        tempvar yhat_n XB_n
        predict `yhat_n', xb
        gen `XB_n' = `yhat_n' - `cd'*`d' - `lam2'*_b[`lam2']
        gen `coveu' = _b[`lam2']

        tempname imenor_mean emenor_mean
        egen `imenor_mean' = mean(`cd'*`prob')
        scalar imenor = `imenor_mean'
        egen `emenor_mean' = mean(`adnadwW'*`dens'*(`XB_n'-`coveu'*`wpro'))
        scalar emenor = `emenor_mean'
        scalar totnor = imenor + emenor

        /* Quadratic model */
        gen `lam3' = 1-`lam2'*`dwpro'
        reg `y' `dx' `lam2' `lam3' if `q'==1
        gen `cd2' = _b[`d']
        foreach var of local outcome_int {
            replace `cd2' = `cd2' + _b[`dvar_map_out_`var''] * `var'
        }
        tempvar yhat_q XB_q
        predict `yhat_q', xb
        gen `XB_q' = `yhat_q' - `cd2'*`d' - `lam2'*_b[`lam2'] - `lam3'*_b[`lam3']
        gen `g1' = _b[`lam2']
        gen `g2' = _b[`lam3']

        tempname imequa_mean emequa_mean
        egen `imequa_mean' = mean(`cd2'*`prob')
        scalar imequa = `imequa_mean'
        egen `emequa_mean' = mean(`adnadwW'*`dens'*(`XB_q'-`g1'*`wpro'+`g2'*`wpro'^2))
        scalar emequa = `emequa_mean'
        scalar totqua = imequa + emequa

        /* Cubic model */
        gen `lam4' = `lam2' * (2 + `dwpro'^2)
        reg `y' `dx' `lam2' `lam3' `lam4' if `q'==1
        gen `cd3' = _b[`d']
        foreach var of local outcome_int {
            replace `cd3' = `cd3' + _b[`dvar_map_out_`var''] * `var'
        }

        tempvar yhat_c XB_c
        predict `yhat_c', xb
        gen `XB_c' = `yhat_c' - `cd3'*`d' - `lam2'*_b[`lam2'] - `lam3'*_b[`lam3'] - `lam4'*_b[`lam4']

        replace `g1' = _b[`lam2']
        replace `g2' = _b[`lam3']
        gen     `g3' = _b[`lam4']

        tempname imecub_mean emecub_mean
        egen `imecub_mean' = mean(`cd3'*`prob')
        scalar imecub = `imecub_mean'
        egen `emecub_mean' = ///
            mean(`adnadwW'*`dens'*(`XB_c' - `g1'*`wpro' + `g2'*`wpro'^2 - `g3'*`wpro'^3))
        scalar emecub = `emecub_mean'
        scalar totcub = imecub + emecub
    }
    
    // Return results
    foreach scalar in tob ///
        imeind emeind totind ///
        imenor emenor totnor ///
        imequa emequa totqua ///
        imecub emecub totcub {
        return scalar `scalar' = scalar(`scalar')
    }
    
    matrix vif = (tob, ///
                  imeind, emeind, totind, ///
                  imenor, emenor, totnor, ///
                  imequa, emequa, totqua, ///
                  imecub, emecub, totcub)
    matrix colnames vif = tob ///
        imeind emeind totind ///
        imenor emenor totnor ///
        imequa emequa totqua ///
        imecub emecub totcub
    return matrix vif = vif
end
