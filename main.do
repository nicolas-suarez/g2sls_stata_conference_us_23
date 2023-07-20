version 18
clear all
set more off, permanently
cd "${dropbox}/peer effects/stata_conference"


*loading adjacency matrix G (4434 by 4434)
mata: mata matuse "..\adjacency.mmat", replace

*loading and preparing data
use "..\base_final_students.dta", clear
lab var gpa_first "First semester GPA"
lab var adm_score "Admission score"
rename sipee aff_action 
lab var aff_action "Affirmative action"
lab var female "Female"
gen major_econ=(cod_mencion=="ME")
lab var major_econ "Major in Economics"
gen major_buss=(cod_mencion=="MA")
lab var major_buss "Major in Business"

*first snippet of code
sjlog using "logs/log1", replace
describe gpa_first adm_score aff_action female major*
list gpa_first adm_score aff_action female major* in 1/5
sjlog close, replace

*base IV model
sjlog using "logs/log2", replace
g2sls gpa_first female aff_action adm_score, row adj(G)
sjlog close, replace
est sto m4

*IV model with fixed effects
sjlog using "logs/log3", replace
g2sls gpa_first female aff_action adm_score, row adj(G) fixed
sjlog close, replace
est sto m5

*OLS model
sjlog using "logs/log4", replace
g2sls gpa_first female aff_action adm_score, row adj(G) ols
sjlog close, replace
est sto m1

*direct effects
sjlog using "logs/log5", replace
g2sls gpa_first female aff_action adm_score, row adj(G) directvariables(major_*)
sjlog close, replace

*other models for table
g2sls gpa_first female aff_action adm_score, row adj(G) ols fixed
est sto m2

g2sls gpa_first female aff_action adm_score, row adj(G) ols fixed dir(major_*)
est sto m3

g2sls gpa_first female aff_action adm_score, row adj(G) fixed directvariables(major_*)
est sto m6

*exporting tables
estout m1 m2 m3 m4 m5 m6 using "tex/results.tex", ///
replace style(tex) collabels(, none) mlabels(, none) ///
stats(N , fmt(%9.0gc) labels("\midrule \textbf{Observations}" ) ) noomitted ///
cells(b(star fmt(%9.4f)) se(par))  starl(* 0.10 ** 0.05 *** 0.01) substitute("(.)" "" )  ///
order(gpa_first_p female_p aff_action_p adm_score_p female aff_action adm_score major*) ///
varlabel(_cons "Constant" gpa_first_p "GPA of peers" female "Female" aff_action "Affirmative Action program" /// 
	adm_score "Admission score" female_p "Share of female peers" aff_action_p "Share of peers in Aff. Action program" adm_score_p "Adm. Score of peers") ///
 label


