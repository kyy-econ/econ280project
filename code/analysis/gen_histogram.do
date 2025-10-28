clear 
cls
cd "C:\Users\Qone9\Desktop\econ280project\code\analysis"

use "C:\Users\Qone9\Desktop\econ280project\data\rawdata\MJLEEdata.dta", clear

* 10 equal-width bins
histogram salarios, percent bin(10) ///
    title("Salarios: Histogram with 10 equal-width bins")

* ensure output folder exists, then export
cap mkdir "C:\Users\Qone9\Desktop\econ280project\output"
graph export "C:\Users\Qone9\Desktop\econ280project\output\salarios_hist_equalwidth.png", replace
graph save  "C:\Users\Qone9\Desktop\econ280project\output\salarios_hist_equalwidth.gph", replace