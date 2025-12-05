program define MJleeIMEEMEs2, rclass
    version 14.0
    syntax varlist(min=3 max=3) ///
        [, ///
        SELection_int(string) ///
        SELection_cov(string) ///
        OUTcome_int(string)   ///
        OUTcome_cov(string)   ///
        REPetitions(integer 50) ///
        SEED(integer -1)]

    display _n(2) as result ///
        "Extensive and Intensive Margin Effects in Sample Selection Models: Racial Effects on Wages, Myoung-jae Lee (2017)" _n ///
        as result "— code by Kyoowon Yang, originally created on 2024/10/4, last revised 2024/10/15. version 1.0.3 (with cubic)." _n(3)
    
    * Parse the three required variables
    tokenize `varlist'
    local labhr    "`1'"
    local emptrain "`2'"
    local worked   "`3'"
    
    * Set seed if specified
    if `seed' >= 0 {
        set seed `seed'
    }
    
    tempfile temp_data
    tempvar temp_est
    
    display _n(2) as text ///
        "Covariates and selection correction terms for probit and tobit estimators and LSEs" _n
    
    *------------------------------------------------------*
    * 1. Raw regressions (probit/tobit/OLS/Heckit/quad/cub)*
    *------------------------------------------------------*
    qui {
        * Generate basic variables
        gen d = `emptrain'
        gen q = `worked'
        gen y = `labhr'
        
        * Generate interaction variables with proper names
        local sel_int_d_f ""
        foreach var of local selection_int {
            gen d`var' = `var' * d
            local sel_int_d_f "`sel_int_d_f' d`var'"
        }
        
        local out_int_d_f ""
        foreach var of local outcome_int {
            if strpos("`sel_int_d_f'", "d`var'") == 0 {
                gen d`var' = `var' * d
            }
            local out_int_d_f "`out_int_d_f' d`var'"
        }
        
        * Create regressor lists
        local dw "d `selection_cov' `sel_int_d_f'"
        local dx "d `outcome_cov' `out_int_d_f'"
        
        * Create ordered variable list for display - put 'd' variables first
        local all_vars "`dw' `dx'"
        local myorder "d"
        
        * First add all variables starting with 'd'
        foreach var of local all_vars {
            if substr("`var'", 1, 1) == "d" & "`var'" != "d" {
                local myorder "`myorder' `var'"
            }
        }
        
        * Then add remaining variables
        foreach var of local all_vars {
            if substr("`var'", 1, 1) != "d" {
                local myorder "`myorder' `var'"
            }
        }
        
        *-----------------------------*
        * Run models for display only *
        *-----------------------------*
        probit q `dw'
        estimates store probit1
        
        * Probit index for IMR / control functions
        predict wpro_2, xb
        gen dwpro = wpro_2
        
        * Tobit
        tobit y `dx', ll(0)
        estimates store tobit1
        
        * Independence (OLS on selected)
        reg y `dx' if q==1
        estimates store indp1
        
        * Normality (Heckit)
        gen lam2 = normalden(dwpro) / normal(dwpro)
        reg y `dx' lam2 if q==1
        estimates store norm1
        
        * Quadratic control function
        gen lam3 = 1 - lam2*dwpro
        reg y `dx' lam2 lam3 if q==1
        estimates store quad1
        
        * Cubic control function: lam4 = lam2 * (2 + dwpro^2)
        gen lam4 = lam2 * (2 + dwpro^2)
        reg y `dx' lam2 lam3 lam4 if q==1
        estimates store cubic1
    }
    
    * Display regression results with corrected ordering
    qui estout probit1 tobit1 indp1 norm1 quad1 cubic1, ///
        cells("b se p") order(`myorder')
    qui return list
    noi matlist r(coefs)
    
    * Drop generated variables before proceeding with bootstrap
    cap drop d
    cap drop q
    cap drop y
    cap drop dwpro
    cap drop wpro_2
    cap drop lam2
    cap drop lam3
    cap drop lam4
    foreach var of local selection_int {
        cap drop d`var'
    }
    foreach var of local outcome_int {
        cap drop d`var'
    }
    
    display _n(2) as text ///
        "Treatment effects (and non-parametric bootstrap 95% confidence intervals)" _n
        
    *---------------------------------------------------------*
    * 2. Call MJleeIMEEME_BASE2 once to get full-sample values*
    *---------------------------------------------------------*
    qui {
        MJleeIMEEME_BASE2 `labhr' `emptrain' `worked', ///
            selection_int(`selection_int') ///
            selection_cov(`selection_cov') ///
            outcome_int(`outcome_int') ///
            outcome_cov(`outcome_cov')
    }
    
    * Store full sample estimates (13 scalars, including cubic)
    matrix define full_est = ///
        (r(tob), ///
         r(imeind), r(emeind), r(totind), ///
         r(imenor), r(emenor), r(totnor), ///
         r(imequa), r(emequa), r(totqua), ///
         r(imecub), r(emecub), r(totcub))
    
    *--------------------------------------------*
    * 3. Manual nonparametric bootstrap of TE's  *
    *--------------------------------------------*
    
    * Create matrix to store bootstrap results (repetitions × 13)
    matrix define boot_results = J(`repetitions', 13, .)
    
    * Start bootstrap
    display "Starting bootstrap with `repetitions' replications..."
    
    * Run bootstrap manually
    forvalues i = 1/`repetitions' {
        preserve
        bsample
        
        qui {
            MJleeIMEEME_BASE2 `labhr' `emptrain' `worked', ///
                selection_int(`selection_int') ///
                selection_cov(`selection_cov') ///
                outcome_int(`outcome_int') ///
                outcome_cov(`outcome_cov')
            
            matrix boot_results[`i',1] = ///
                (r(tob), ///
                 r(imeind), r(emeind), r(totind), ///
                 r(imenor), r(emenor), r(totnor), ///
                 r(imequa), r(emequa), r(totqua), ///
                 r(imecub), r(emecub), r(totcub))
        }
        restore
        
        display "." _continue
        if mod(`i', 10) == 0 {
            display "`i'"
        }
    }
    display
    
    *-----------------------------------*
    * 4. Compute bootstrap SEs (1×13)   *
    *-----------------------------------*
    matrix acc_results = J(1, 13, 0)
    matrix acc_squared = J(1, 13, 0)
    
    forvalues i = 1/`repetitions' {
        matrix this_row = boot_results[`i',1...]
        matrix acc_results = acc_results + this_row
        
        matrix this_squared = this_row
        matrix this_squared = hadamard(this_squared, this_squared)
        matrix acc_squared = acc_squared + this_squared
    }
    
    matrix means         = acc_results / `repetitions'
    matrix means_squared = hadamard(means, means)
    matrix squared_means = acc_squared / `repetitions'
    matrix var           = squared_means - means_squared
    
    * Standard errors: sqrt of variances
    matrix se = J(1, 13, 0)
    forvalues i = 1/13 {
        matrix se[1,`i'] = sqrt(var[1,`i'])
    }
    
    *---------------------------------------------*
    * 5. Percentile-based 95% confidence intervals*
    *---------------------------------------------*
    matrix ci_lower = J(1, 13, 0)
    matrix ci_upper = J(1, 13, 0)
    
    qui {
        preserve
        clear
        
        * Create dataset from bootstrap results
        set obs `repetitions'
        
        forvalues j = 1/13 {
            gen v`j' = .
            forvalues i = 1/`repetitions' {
                replace v`j' = boot_results[`i',`j'] in `i'
            }
            
            _pctile v`j', p(2.5 97.5)
            matrix ci_lower[1,`j'] = r(r1)
            matrix ci_upper[1,`j'] = r(r2)
        }
        restore
    }
    
    *-----------------------------------------------------*
    * 6. Assemble final results table (13 × 4) and return *
    *-----------------------------------------------------*
    matrix results = (full_est' , se', ci_lower', ci_upper')
    matrix colnames results = "Estimate" "Std.Err." "95% CI Lower" "95% CI Upper"
    matrix rownames results = ///
        "Tobit"   ///
        "IME_ind" "EME_ind" "TE_ind" ///
        "IME_nor" "EME_nor" "TE_nor" ///
        "IME_qua" "EME_qua" "TE_qua" ///
        "IME_cub" "EME_cub" "TE_cub"
    
    * Display final matrix
    display _n ///
        "Estimates, Standard Errors, and 95% Percentile Bootstrap Confidence Intervals (`repetitions' replications)" _n
    matrix list results, format(%9.4f)
    
    * Return results
    return matrix results           = results
    return matrix full_estimates    = full_est
    return matrix bootstrap_results = boot_results
    return matrix standard_errors   = se
    return matrix ci_lower          = ci_lower
    return matrix ci_upper          = ci_upper
end
