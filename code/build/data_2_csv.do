*******************************************************
* 01_export_raw_to_csv.do
* 
* Purpose:
*   - Load the original wage sample used in Lee (2017)
*   - Export it as a CSV file for transparency and portability
*
* Note:
*   This paper is an econometrics *theory* paper. The focus of this
*   replication project is not on data cleaning but on implementing
*   and extending the decomposition method.
*
*   Therefore:
*   - We do NOT perform any data cleaning or transformation here.
*   - The CSV file is simply a direct export of the original Stata data.
*
*   The primary purpose of this repository is to provide a Stata command
*   (the MJleeIMEEME* ado files) that can be applied to *general*
*   datasets to obtain the extensive and intensive margin effects
*   reported in the paper.
*******************************************************

clear all
set more off

use "MJLEEdata.dta", clear
export delimited using "MJLEEdata.csv", replace
