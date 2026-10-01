Rcpp::NumericVector compute_logC(const unsigned int& n, const double& scale, const double& location){

	if(!( (scale<0) & (location<=0) ) )
		throw std::runtime_error("Error in my_logC. The recursive formula for the absolute values of the C numbers can be you used if the scale is strictly negative and location in non positive");

	const double& s = -scale; //s is strictly positive
	const double& r = -location; //r is non-negative


	Rcpp::NumericVector LogC_old(n+1, 0.0);

	if(n == 0)
		return LogC_old; // nothing to do in this case

	// Compute the first row
	LogC_old[0] = log_raising_factorial(1,r);
	LogC_old[1] = std::log(s);

	//Rcpp::NumericVector LogC_update(LogC_old);
	Rcpp::NumericVector LogC_update(n+1, 0.0);
	double coef(0.0);
	for(std::size_t nn = 2; nn <= n; ++nn ){ //for each row
		LogC_update[0] = log_raising_factorial(nn,r);
		for(std::size_t k = 1; k < nn; ++k){ //for each column but the first and the last one
			coef = s*k + r + nn - 1;
			LogC_update[k] = std::log(coef) + LogC_old[k] + std::log( 1 + s/coef*std::exp( LogC_old[k-1] - LogC_old[k] ) );
		}
		LogC_update[nn] = nn*std::log(s); //update last element

		//std::copy(LogC_update.begin(),LogC_update.end(),LogC_old.begin());
		std::swap(LogC_old, LogC_update); //avoid copy, LogC_update has to be overwritten but in this way no copy is performed to update LogC_old.
	}
	return (LogC_old);
}