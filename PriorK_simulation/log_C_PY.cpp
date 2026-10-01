#include <Rcpp.h>
#include <cmath>
#include <limits>

double log_add_py(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) return b;
    if (b == -std::numeric_limits<double>::infinity()) return a;
    double hi = std::max(a, b);
    return hi + std::log1p(std::exp(std::min(a, b) - hi));
}

// S(n,k;sigma) = C_mathcal(n,k;sigma)/sigma^k.
// S(n,k) = (n-1-sigma*k) S(n-1,k) + S(n-1,k-1).
// sigma=0: limite esatto, Stirling non segnati di prima specie.
// Output: log S(n,k), k=0,...,n (indice R k+1).
// [[Rcpp::export]]
Rcpp::NumericVector log_C_scaled_PY(int n, double sigma) {
    if (n < 0 || !std::isfinite(sigma) || sigma < 0 || sigma >= 1)
        Rcpp::stop("Richiesti n>=0 e 0<=sigma<1.");
    const double neginf = -std::numeric_limits<double>::infinity();
    Rcpp::NumericVector old(n+1, neginf);
    old[0] = 0.0;
    for (int m = 1; m <= n; ++m) {
        Rcpp::NumericVector next(n+1, neginf);
        for (int k = 1; k <= m; ++k) {
            double same = neginf;
            if (k < m) same = std::log((m-1)-sigma*k) + old[k];
            next[k] = log_add_py(same, old[k-1]);
        }
        old = next;
    }
    return old;
}

// C_mathcal di Lijoi, Mena e Pruenster (2007), appendice p.782.
// C(n,k) = (n-1-sigma*k) C(n-1,k) + sigma*C(n-1,k-1).
// Nessun fattore di segno aggiuntivo nella PMF.
// [[Rcpp::export]]
Rcpp::NumericVector log_C_PY(int n, double sigma) {
    Rcpp::NumericVector ans = log_C_scaled_PY(n, sigma);
    for (int k = 1; k <= n; ++k)
        ans[k] = sigma == 0 ? -std::numeric_limits<double>::infinity()
                            : ans[k] + k * std::log(sigma);
    return ans;
}
