#include <Rcpp.h>
#include <cmath>
#include <limits>

// Dipendenza mancante nel file originale: log del fattoriale crescente.
double log_raising_factorial(const unsigned int& n, const double& a)
{

    if(n==0)
        return 0.0;
    if(a<0)
        throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale");
    else if(a==0.0){
        return -std::numeric_limits<double>::infinity();
    }
    else{

        double val_max{std::log(a+n-1)};
        double res{1.0};
        if (n==1)
            return val_max;
        for(std::size_t i = 0; i <= n-2; ++i){
            res += std::log(a + (double)i) / val_max;
        }
        return val_max*res;

    }
}

// Include direttamente l'implementazione fornita, senza modificarla.
#define ORIGINAL_LOGC_FILE "../compute_logC.cpp"
#include ORIGINAL_LOGC_FILE

// [[Rcpp::export]]
Rcpp::NumericVector original_logC(unsigned int n, double scale, double location = 0) {
    return compute_logC(n, scale, location);
}
