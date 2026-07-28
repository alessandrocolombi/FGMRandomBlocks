##' get group indexes
#'
#' @param rho partition in che compact version ef c(1,3,3)
#'
#' @return a vector whose components are
#' the increasing indexes of the elements where each group starts
#' e.g. if rho=c(1,3,3), the output will be c(1,4,6),
#' meaning that 1 is the index of the first element of the first group,
#' 4 is the index of the first element of the second group,
#' 6 is the index of the first element of the third group
#'
#' @export
#'
#' @examples
get_group_indexes = function(rho){
    return(cumsum(rho))
}

# =============================================================================
# =           FUNCIONS FOR MOVING BETWEEN PARTITION REPRESENTATIONS           =
# =============================================================================

# Vecchia versione del Poli. Buggata se rho = p

# rho_to_r = function(rho){
#     group_indexes = get_group_indexes(rho)
#     group_indexes = group_indexes[1:(length(group_indexes)-1)]
#     r <- numeric(sum(rho)-1)
#     r[group_indexes] = 1
#     return(r)
# }

rho_to_r = function(rho){
  if(length(rho) == 1)
    return( rep(0,rho-1) )
  
  group_indexes = get_group_indexes(rho)
  group_indexes = group_indexes[1:(length(group_indexes)-1)]
  r <- numeric(sum(rho)-1)
  r[group_indexes] = 1
  return(r)
}

z_to_rho = function(z){
    return(as.vector(table(z)))
}

z_to_r = function(z){
    return(diff(z))
}

rho_to_z = function(rho){
    z = numeric()
    
    for(i in seq_along(rho)){
        z = c(z, rep(i, rho[i]))
    }
    return(z)
}


# Luciano Version
r_to_rho <- function(rvec) {
  cp <- which(rvec == 1)         # changepoint positions
  # if (length(cp) == 0) return(integer(0))   # no changepoints then empty rho
  rho <- diff(c(0, cp))          # segment lengths
  return(rho)
}

r_to_z = function(r){
    if(length(r) == 0)
        return(integer(0))
    
    r = as.integer(r)
    
    if(any(is.na(r)) || any(!r %in% c(0L, 1L)))
        stop("r must be a binary changepoint vector")
    
    if(r[length(r)] != 1L)
        r = c(r, 1L)
    
    rho_to_z(r_to_rho(r))
}



#' Compute the rising factorial (also called Pochhammer symbol)
#'
#' Peturns the value of Pochhammer's symbol calculated as
#' \deqn{(x)_n = x (x+1) \cdots (x+n-1)}
#'
#'
#' @param x numeric value for the argument of the symbol
#' @param n integer value for the number of terms in the symbol
#' @param log boolean value, if TRUE the rising factorial is returned in log. Default is TRUE.
#'
#' @return the rising factorial of x with n terms.
#' @export
#'
#' @examples
# lpochhammer <- function(x, n, log = TRUE) {
#     if (n < 0)
#         stop("The pochhammer operator doesn't allow n < 0")
#     if (x < 0)
#         stop("The pochhammer operator doesn't allow x < 0")
#     if (x == 0){
#         if(log)
#             return(NA)
#         else
#             return(0)
#     }
# 
#     num_vec <- numeric(n)
#     
#     if (n != 0) {
#         num_vec[1] = x
#         for (i in 1:(n - 1)) {
#             num_vec[i + 1] = (x + i)
#         }
#     }
#     
#     if (n == 0) {
#         num_vec[1] = 1
#     }
#     
#     if (log)
#         return(sum(log(num_vec)))
#     else
#         return(prod(num_vec))
# }

# Luciano version
lpochhammer <- function(x, n, log = TRUE) {
  if (n < 0) stop("n must be >= 0")
  
  if (log) {
    return(lgamma(x + n) - lgamma(x))
  } else {
    return(exp(lgamma(x + n) - lgamma(x)))
  }
}



#' Absolute value of Stirling number of the first kind (adapted)
#' Computes an adapted version of the Stirling number of the first kind
#'
#' The Stirling number represents the number of ways that we can arrange
#' k objects around indistinguishable circles of length j
#'
#' NOTE: Requires gmp package.
#'
#' @param k First parameter - indicates the overall number of objects
#' @param j Second parameter - indicates the length of the circles (see above)
#'
#' @return a positive scalar indicating the adapted version of the Stirling number of the first kind
#' (i.e. the "unacceptable" values are turned to zeroes)
#' @export
#'
#' @examples
abs_stirling_number_1st <- function(k,j){
    if (j == 0 && k == 0) {
        return(1)
    }
    
    if (k < 0) {
        stop("In computing the Stirling number of the first kind, k must be greater or equal than 0.")
    }
    if (j <= 0 || j > k) {
        abs_stir_num_first_kind = 0
    }
    else{
        abs_stir_num_first_kind = (as.numeric(abs(gmp::Stirling1(k, j))))
    }
    return(abs_stir_num_first_kind)
}


#' Shifted gamma function
#'
#' Computes a shifted gamma, given the parameters and the shift.
#' Z ~ shiftedGamma(alpha,beta,mu) is and only if Z-mu ~ Gamma(alpha,beta)
#'
#' @param alpha first parameter of the Gamma
#' @param beta second parameter of the Gamma
#' @param mu shift parameter
#'
#' @return
#' @export
#'
#' @examples
shifted_gamma <- function(alpha, beta, mu, n = 1) {
    rgamma(n, alpha, beta) + mu
}



#' Bayesian FDR Analysis
#'
#' \loadmathjax Given the plinks matrix, this utility computes the False Discovery Rate Index, forcing the false discovery rates to be less than \code{min_rate.}
#' @param plinks matrix containing the posterior inclusion probability for each link. It has to be upper triangular. Its dimension depends on the type of graph it represents.
#' It is indeed possible to pass a \mjseqn{p \times p} matrix, or a \mjseqn{n\_groups \times n\_groups}.
#' @param tol sequence of tolerances to be tested trying to select a graph truncating \code{plinks} at that value.
#' @param min_rate fix false discoveries to remain under this selected threshold.
#' @param diag boolean, if the diagonal of \code{plinks} has to be included in the computations. Set \code{FALSE} if the graph is in complete form, set \code{TRUE} for block graphs.
#'
#' @return a list of two elements: best_threshold, the best value of tol according to this analysis.
#' best_truncated_graph, the proposed posterior graph according to the analysis.
#' @export
BFDR_selection = function (plinks, tol = seq(0.1, 1, by = 0.025), min_rate = 0.05, diag = FALSE)
{
    if(dim(plinks)[1] != dim(plinks)[2])
        stop("plinks matrix must to be squared")
    p = dim(plinks)[1]
    plinks_vet = plinks[upper.tri(plinks, diag = diag)]
    if (any(tol > max(plinks_vet)))
        tol <- tol[-which(tol > max(plinks_vet))]
    if (is.null(tol))
        stop("No feasible tolerances")
    FDR = rep(0, length(tol))
    for (i in 1:length(tol)) {
        tolerance <- tol[i]
        above_tr = plinks_vet[plinks_vet >= tolerance]
        FDR[i] = sum(1 - above_tr)/length(above_tr)
    }
    if(FDR[1] < min_rate) {
        best_soglia_fdr = tol[1]
    }
    else for (i in 2:length(FDR)) {
        if (FDR[i] < min_rate)
            (break)()
    }
    best_soglia_fdr = tol[i]
    best_graph_fdr = matrix(0, p, p)
    best_graph_fdr[plinks >= best_soglia_fdr] = 1
    result = list(
        "best_treshold"=best_soglia_fdr,
        "best_truncated_graph"=best_graph_fdr
    )
    return(result)
}


ACheatmap = function(Mat, 
                     center_value = 0.5, col.upper = "#6D0026", col.center = "#FFBFAA", col.lower = "white",
                     col.n_breaks = 59, use_x11_device = TRUE, remove_diag = FALSE, main = " ", x_label = " ",
                     y_label = " ", horizontal=TRUE )
{
  
  #library(fields)
  #Check for NA
  if(any(is.na(Mat))){
    cat('\n NA values have been removed from Matrix  \n')
  }
  
  #Create color palette
  if(col.n_breaks %% 2 == 0){
    warning( 'col.n_breaks is even but it has to be odd. Adding 1' )
    col.n_breaks = col.n_breaks + 1
  }
  colorTable = fields::designer.colors(col.n_breaks, c( col.lower, col.center, col.upper) )
  col_length = (col.n_breaks + 1) / 2
  
  #Check Matrix
  if(remove_diag){
    diag(Mat) = NA
  }
  min_val = min(Mat,na.rm=T)
  max_val = max(Mat,na.rm=T)
  p_row = dim(Mat)[1]
  p_col = dim(Mat)[2]
  
  
  #Plot
  if(!is.null(center_value)){
    
    if(!(min_val < center_value & max_val > center_value)){
      stop('\n The lowest value has to be smaller than center_value and the highest value has to be larger. \n')
    }
    brks = c(seq( min_val, center_value-0.0001,l=col_length), seq( center_value+0.0001, max_val, l=col_length))
    
  }else{
    brks = seq( min_val, max_val, l=2*col_length)
  }
  
  colnames(Mat) = 1:p_col
  rownames(Mat) = 1:p_row
  
  if(use_x11_device){
    x11()
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1),mgp=c(3,1,0))
  if(horizontal){
    
    fields::image.plot(Mat, axes=F, horizontal=T, main = main,
                       col=colorTable,breaks=brks,xlab=x_label,ylab=y_label)
  }else{
    
    fields::image.plot(Mat, axes=F, horizontal=FALSE, main = main,
                       col=colorTable,breaks=brks,xlab=x_label,ylab=y_label)
  }
  box()
  
}





ACheatmap_nolegend = function(Mat, 
                     center_value = 0.5, col.upper = "#6D0026", col.center = "#FFBFAA", col.lower = "white",
                     col.n_breaks = 59, use_x11_device = TRUE, remove_diag = FALSE, main = " ", x_label = " ",
                     y_label = " ")
{
  
  #library(fields)
  #Check for NA
  if(any(is.na(Mat))){
    cat('\n NA values have been removed from Matrix  \n')
  }
  
  #Create color palette
  if(col.n_breaks %% 2 == 0){
    warning( 'col.n_breaks is even but it has to be odd. Adding 1' )
    col.n_breaks = col.n_breaks + 1
  }
  colorTable = fields::designer.colors(col.n_breaks, c( col.lower, col.center, col.upper) )
  col_length = (col.n_breaks + 1) / 2
  
  #Check Matrix
  if(remove_diag){
    diag(Mat) = NA
  }
  min_val = min(Mat,na.rm=T)
  max_val = max(Mat,na.rm=T)
  p_row = dim(Mat)[1]
  p_col = dim(Mat)[2]
  
  
  #Plot
  if(!is.null(center_value)){
    
    if(!(min_val < center_value & max_val > center_value)){
      stop('\n The lowest value has to be smaller than center_value and the highest value has to be larger. \n')
    }
    brks = c(seq( min_val, center_value-0.0001,l=col_length), seq( center_value+0.0001, max_val, l=col_length))
    
  }else{
    brks = seq( min_val, max_val, l=2*col_length)
  }
  
  colnames(Mat) = 1:p_col
  rownames(Mat) = 1:p_row
  
  if(use_x11_device){
    x11()
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1),mgp=c(3,1,0))
  image(
      x = seq(0, 1, length.out = ncol(Mat)),
      y = seq(0, 1, length.out = nrow(Mat)),
      z = Mat,
      axes = FALSE,
      col = colorTable,
      breaks = brks,
      main = main,
      xlab = x_label,
      ylab = y_label
  )
  box()
  
}






# keep_last <- function(chains, n = 5000) {
#   lapply(chains, function(x) {
#     len <- length(x)
#     if (len > n) {
#       x[(len - n + 1):len]
#     } else {
#       x
#     }
#   })
# }

keep_last <- function(chains, n = 5000) {
  lapply(chains, function(x) {
    
    # ---- 3D arrays (e.g. K, graph_samples, Beta)
    if (is.array(x) && length(dim(x)) == 3) {
      d <- dim(x)
      if (d[3] > n) {
        return(x[, , (d[3] - n + 1):d[3], drop = FALSE])
      } else {
        return(x)
      }
    }
    
    # ---- matrices (e.g. mu)
    if (is.matrix(x)) {
      if (nrow(x) > n) {
        return(x[(nrow(x) - n + 1):nrow(x), , drop = FALSE])
      } else {
        return(x)
      }
    }
    
    # ---- vectors
    if (is.vector(x) && !is.list(x)) {
      if (length(x) > n) {
        return(x[(length(x) - n + 1):length(x)])
      } else {
        return(x)
      }
    }
    
    # ---- lists (rho, gamma)
    if (is.list(x)) {
      if (length(x) > n) {
        return(x[(length(x) - n + 1):length(x)])
      } else {
        return(x)
      }
    }
    
    # fallback
    return(x)
  })
}



# keep_thin <- function(chains, burn_in = 0, thin = 10) {
#   lapply(chains, function(x) {
#     niter <- length(x)
#     x[seq(burn_in + 1, niter, by = thin)]
#   })
# }
