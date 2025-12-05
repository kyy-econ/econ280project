cls
clear all
set more off

use MJLEEdata, clear

export delimited using "MJLEEdata.csv", replace

MJleeIMEEMEs2 salarios etnia wagesmpl, ///
    selection_int(age educa voca veterano hisppor calif newme texas marital) ///
    selection_cov(age age2_100 age_edu educa voca veterano hisppor calif newme texas marital) ///
    outcome_int(age educa voca veterano dindphis doccphis hisppor calif newme texas) ///
    outcome_cov(age age2_100 age_edu educa voca veterano dindphis doccphis hisppor calif newme texas) ///
    repetitions(200) seed(12345)


matrix R = r(results)

esttab matrix(R) using "lee_table.tex", replace ///
     b(%9.4f) booktabs ///
     collabels("Estimate" "Std. Err." "95% CI lower" "95% CI upper") ///
     nonumber noobs compress ///
     title("Extensive and Intensive Margin Effects (Lee, 2017)")