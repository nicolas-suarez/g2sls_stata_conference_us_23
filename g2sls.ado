*! version 1.01 20jun2023

capture program drop g2sls
program define g2sls, eclass
	version 18
	syntax varlist [if] [in] , ADJacency(name) [ROW FIXED OLS DIRectvariables(varlist) level(integer 95)] 
	/* implements the generalized 2SLS model of Bramoulle et al (2009), without fixed effects.
	The model includes a constant, then the endogenous effect, effects of the independent variables, and then the exogenous effects.
	Adapted from https://rpubs.com/Nicolas_Gallo/549370

	INPUTS:
	varlist = exogenous variables
	dep= dependent variable
	adjacency=name of the adjacency matrix in Mata
	row=optional, row normalizes the adjacency matrix
	fixed=optional, estimates model with cohort level fixed effects
	ols=optional, OLS results. Don't use together with FIXED
	dirvar=optional. Variables to use in the regression that won't have peer effects
	level=optional, level for the reported CIs

	OUTPUT:
	displays the coefficients and standard errors in a table. Stores eclass results.
	*/
	*separating dependent from independent variables
	gettoken dep varlist: varlist
	preserve
		quietly{

			// *returning error if OLS and FIXED are used together
			// if "`ols'"!="" & "`fixed'"!="" { 
			// 	noisily display as error "options OLS and FIXED may not be combined"
			// 	exit 184
			// }

			// *returning error if VARDIR and OLS are not used together
			// if "`ols'"==""  & "`directvariables'"!="" { 
			// 	noisily display as error "option VARDIR has to be used with OLS option"
			// 	exit 184
			// }

			*marking our sample
			marksample touse
			markout `touse' `dep' `varlist' `directvariables'
			*moving data as matrices
			mata X=st_data(.,"`varlist'")
			mata y=st_data(.,"`dep'")
			
			mata muestra=st_data(.,st_local("touse"))
			
			if "`directvariables'"!="" mata dirvar=st_data(.,"`directvariables'")
			
			*dropping missing values from data matrices
			mata X=select(X,muestra)
			mata y=select(y,muestra)
			if "`directvariables'"!="" mata dirvar=select(dirvar,muestra)
			*dropping missing values from G matrix (eliminating the rows and columns with missing values, so the matrix are comformable)
			mata peer_iv_adj=select(`adjacency',muestra)
			mata peer_iv_adj=select(peer_iv_adj,muestra')
			*row normalizing G if needed
			if "`row'"!="" mata peer_iv_adj=peer_iv_adj:/editvalue(rowsum(peer_iv_adj),0,1)

			*generating identity matrix
			mata Id=I(rows(peer_iv_adj))
			
			*OLS results
			if "`ols'"!="" {
				*with fixed effects
				if "`fixed'"!="" {
					if "`directvariables'"!="" mata X_1 = ( (Id-peer_iv_adj)*peer_iv_adj*y, (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X, (Id-peer_iv_adj)*dirvar)
					else mata X_1 = ( (Id-peer_iv_adj)*peer_iv_adj*y, (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X )
					mata theta= invsym(quadcross(X_1, X_1))*quadcross(X_1, (Id-peer_iv_adj)*y)
					mata e= (Id-peer_iv_adj)*y - X_1*theta
					mata V = (quadsum(e:^2)/(rows(X_1)-cols(X_1)))*invsym(quadcross(X_1, X_1))	
				}
				else{
					if "`directvariables'"!="" mata X_1 = ( J(rows(X),1,1), peer_iv_adj*y, X, peer_iv_adj*X, dirvar)
					else mata X_1 = ( J(rows(X),1,1), peer_iv_adj*y, X, peer_iv_adj*X )
					mata theta= invsym(quadcross(X_1, X_1))*quadcross(X_1, y)
					mata e= y - X_1*theta
					mata V = (quadsum(e:^2)/(rows(X_1)-cols(X_1)))*invsym(quadcross(X_1, X_1))
				}
				
			}
			else {
				*putting matrices together
				*with fixed effects
				if "`fixed'"!="" {
					if "`directvariables'"!="" mata S= (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X, (Id-peer_iv_adj)*peer_iv_adj*peer_iv_adj*X, (Id-peer_iv_adj)*dirvar
					else mata S= (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X, (Id-peer_iv_adj)*peer_iv_adj*peer_iv_adj*X 
					if "`directvariables'"!="" mata X_1= (Id-peer_iv_adj)*peer_iv_adj*y, (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X, (Id-peer_iv_adj)*dirvar
					else mata X_1= (Id-peer_iv_adj)*peer_iv_adj*y, (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X 				
				}
				else{
					if "`directvariables'"!="" mata S= (J(rows(X),1,1), X, peer_iv_adj*X, peer_iv_adj*peer_iv_adj*X, dirvar )
					else mata S= (J(rows(X),1,1), X, peer_iv_adj*X, peer_iv_adj*peer_iv_adj*X )
					if "`directvariables'"!="" mata X_1= (J(rows(X),1,1), peer_iv_adj*y, X, peer_iv_adj*X, dirvar)
					else mata X_1= (J(rows(X),1,1), peer_iv_adj*y, X, peer_iv_adj*X )
				}
				mata P= S*invsym(quadcross(S,S))*S'

				*first 2sls
				if "`fixed'"!="" mata theta_1= invsym(X_1'*P*X_1)*X_1'*P*(Id-peer_iv_adj)*y
				else mata theta_1= invsym(X_1'*P*X_1)*X_1'*P*y
				
				*building instrument
				if "`fixed'"!="" {
					if "`directvariables'"!="" mata Z = peer_iv_adj*luinv(Id-theta_1[1]*peer_iv_adj)*(Id-peer_iv_adj)*(X*theta_1[2::(1+cols(X))] +  peer_iv_adj*X*theta_1[(2+cols(X))::(1+2*cols(X))] + dirvar*theta_1[(2+2*cols(X))::(1+2*cols(X)+cols(dirvar))] ), (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X, (Id-peer_iv_adj)*dirvar
					else mata Z = peer_iv_adj*luinv(Id-theta_1[1]*peer_iv_adj)*(Id-peer_iv_adj)*(X*theta_1[2::(1+cols(X))] +  peer_iv_adj*X*theta_1[(2+cols(X))::(1+2*cols(X))] ), (Id-peer_iv_adj)*X, (Id-peer_iv_adj)*peer_iv_adj*X	
				}
				else{
					if "`directvariables'"!="" mata Z = J(rows(X),1,1), peer_iv_adj*luinv(Id-theta_1[2]*peer_iv_adj)*( theta_1[1]*J(rows(X),1,1) + X*theta_1[3::(2+cols(X))] +  peer_iv_adj*X*theta_1[(3+cols(X))::(2+2*cols(X))] + dirvar*theta_1[(3+2*cols(X))::(2+2*cols(X)+cols(dirvar))]  ), X, peer_iv_adj*X, dirvar
					else mata Z = J(rows(X),1,1), peer_iv_adj*luinv(Id-theta_1[2]*peer_iv_adj)*( theta_1[1]*J(rows(X),1,1) + X*theta_1[3::(2+cols(X))] +  peer_iv_adj*X*theta_1[(3+cols(X))::(2+2*cols(X))] ), X, peer_iv_adj*X
				}
				*
				
				*final IV
				if "`fixed'"!="" mata theta = luinv(quadcross(Z,X_1))*quadcross(Z,(Id-peer_iv_adj)*y)
				else mata theta = luinv(quadcross(Z,X_1))*quadcross(Z,y)

				*resids
				if "`fixed'"!="" {
					if "`directvariables'"!="" mata e= (Id-peer_iv_adj)*y - luinv(Id-theta[1]*peer_iv_adj)*((Id-peer_iv_adj)*X*theta[2::(1+cols(X))] + (Id-peer_iv_adj)*peer_iv_adj*X*theta[(2+cols(X))::(1+2*cols(X))] + (Id-peer_iv_adj)*dirvar*theta[(2+2*cols(X))::(1+2*cols(X)+cols(dirvar))]  )
					else mata e= (Id-peer_iv_adj)*y - luinv(Id-theta[1]*peer_iv_adj)*((Id-peer_iv_adj)*X*theta[2::(1+cols(X))] + (Id-peer_iv_adj)*peer_iv_adj*X*theta[(2+cols(X))::(1+2*cols(X))] )
				}
				else{
					if "`directvariables'"!="" mata e= y - luinv(Id-theta[2]*peer_iv_adj)*( theta[1]*J(rows(X),1,1) + X*theta[3::(2+cols(X))] +  peer_iv_adj*X*theta[(3+cols(X))::(2+2*cols(X))] +dirvar*theta[(3+2*cols(X))::(2+2*cols(X)+cols(dirvar))]  )
					else mata e= y - luinv(Id-theta[2]*peer_iv_adj)*( theta[1]*J(rows(X),1,1) + X*theta[3::(2+cols(X))] +  peer_iv_adj*X*theta[(3+cols(X))::(2+2*cols(X))] )
				}


				*variance
				mata V = luinv(quadcross(Z,X_1))*(Z')*diag(e:^2)*Z*luinv(quadcross(X_1,Z))
			}

			*sending results to Stata
			mata st_matrix("b",theta')
			mata st_matrix("V",V)

			*row and col names for matrices
			local exog_peer //list for names of exogenous effects
			foreach var in `varlist'{
				local exog_peer `exog_peer' `var'_p
			}
			if "`fixed'"!="" {
				local varnames `dep'_p `varlist' `exog_peer' `directvariables'
			}
			else{
				local varnames _cons `dep'_p `varlist' `exog_peer' `directvariables'
			}
			

			*adding col and rownames
			matrix colnames b= `varnames'
			matrix colnames V = `varnames'
			matrix rownames V = `varnames'
		}
		*storing eclass results


		ereturn post b V, depname(`dep') esample(`touse')
		mata st_numscalar("e(N)", rows(peer_iv_adj))
		mata st_numscalar("e(df_r)", rows(X_1)-cols(X_1))
		eret local cmd peer_iv

		di ""
		di _col(50) " Number of obs = " in ye %10.0g e(N)
		if "`fixed'"!="" di "Controlling for component-level fixed effects"
		ereturn display, level(`level')

	restore		
 end
