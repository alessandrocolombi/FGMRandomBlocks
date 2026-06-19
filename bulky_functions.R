library(logr)

#' Main function that updates the partition
#'
#' @param rho The partition in compact form (e.g. rho=c(1,4,5) means that the first group has 1 element, the second has 4 elements and the last has 5 elements).
#' @param alpha_add Fixed probability of choosing an add move or delete move.
#' @param weights_a Vector of size (number of nodes - 1) containing at element j the weights to consider when ADDING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param weights_d Vector of size (number of nodes - 1) containing at element j the weights to consider when DELETING a changepoint between point j and point j+1 (weights are non-normalized probabilities).
#' @param theta_prior Prior parameter as in Martinez and Mena (2014)
#' @param sigma_prior Prior parameter as in Martinez and Mena (2014)
#' @param G Adjacency matrix of the graph
#' @param beta_params Parameters of the Beta
#'
#' @return a new partition in the compact form
#' @export
#'
#' @examples
# rho_current = rho
# G = graph
# update_partition = function(rho_current,
#                             alpha_add,
#                             weights_a,
#                             weights_d,
#                             theta_prior,
#                             sigma_prior,
#                             G,
#                             beta_params) {
#     unifsample = runif(n = 1)
#     choose_add = unifsample < alpha_add
#     
#     # number of groups
#     M = length(rho_current)
#     p = sum(rho_current)
#     
#     # force opposite choice if merging/splitting is not feasible
#     if ((!choose_add && M == 1) || (choose_add && M == p)) {
#         choose_add = !choose_add
#     }
# 
#     proposal_list = proposal_ratio(rho_current, alpha_add, weights_a, weights_d, choose_add)
#     log_proposal_ratioNow = log(proposal_list$ratio)
#     candidate = proposal_list$candidate
#     
#     # compute proposed partition based on candidate and index of the group
#     if (choose_add) {
#         list_output_modify_partition = split_partition(candidate, rho_current)
#     } else {
#         list_output_modify_partition = merge_partition(candidate, rho_current)
#     }
#     rho_proposed        = list_output_modify_partition$new_rho
#     changed_group_index = list_output_modify_partition$changed_group_index
#     
#     #log_print("rho_current", console = FALSE)
#     #log_print(rho_current, console = FALSE)
#     #log_print("rho_proposed", console = FALSE)
#     #log_print(rho_proposed, console = FALSE)
# 
#     log_prior_ratioNow = log_prior_ratio(
#         theta_prior,
#         sigma_prior,
#         rho_current,
#         rho_proposed,
#         choose_add,
#         changed_group_index
#     )
#     
#     log_likelihood_ratioNow = log_likelihood_ratio(
#         alpha_add,
#         weights_a,
#         weights_d,
#         G,
#         rho_current,
#         rho_proposed,
#         choose_add,
#         beta_params$alpha,
#         beta_params$beta,
#         changed_group_index
#     )
#     
#     alpha_accept <- min(1, exp(log_likelihood_ratioNow +
#                                log_prior_ratioNow +
#                                log_proposal_ratioNow))
# 
#     if (runif(n = 1) < alpha_accept) {
#         accepted = TRUE
#         #log_print("Move ACCEPTED", console = FALSE)
#         rho_updated = rho_proposed
#     } else {
#         accepted = FALSE
#         #log_print("Move REJECTED", console = FALSE)
#         rho_updated = rho_current
#     }
#     return(
#         list(
#             "rho_updated" = rho_updated,
#             "accepted" = as.numeric(accepted),
#             "choose_add" = choose_add,
#             "candidate" = candidate,
#             'alpha_accept' = alpha_accept
#         )
#     )
#     
# }
update_partition = function(rho_current,
                            alpha_add,
                            weights_a,
                            weights_d,
                            theta_prior,
                            sigma_prior,
                            G,
                            beta_params) {
  unifsample = runif(n = 1)
  choose_add = unifsample < alpha_add
  
  # number of groups
  M = length(rho_current)
  p = sum(rho_current)
  
  # force opposite choice if merging/splitting is not feasible
  if ((!choose_add && M == 1) || (choose_add && M == p)) {
    choose_add = !choose_add
  }
  
  proposal_list = proposal_ratio(rho_current, alpha_add, weights_a, weights_d, choose_add)
  log_proposal_ratioNow = log(proposal_list$ratio)
  candidate = proposal_list$candidate
  
  # compute proposed partition based on candidate and index of the group
  if (choose_add) {
    list_output_modify_partition = split_partition(candidate, rho_current)
  } else {
    list_output_modify_partition = merge_partition(candidate, rho_current)
  }
  rho_proposed        = list_output_modify_partition$new_rho
  changed_group_index = list_output_modify_partition$changed_group_index
  
  #log_print("rho_current", console = FALSE)
  #log_print(rho_current, console = FALSE)
  #log_print("rho_proposed", console = FALSE)
  #log_print(rho_proposed, console = FALSE)
  
  log_prior_ratioNow = log_prior_ratio(
    theta_prior,
    sigma_prior,
    rho_current,
    rho_proposed,
    choose_add,
    changed_group_index
  )
  
  log_likelihood_ratioNow = log_likelihood_ratio(
    alpha_add,
    weights_a,
    weights_d,
    G,
    rho_current,
    rho_proposed,
    choose_add,
    beta_params$alpha,
    beta_params$beta,
    changed_group_index
  )
  
  alpha_accept <- min(1, exp(log_likelihood_ratioNow +
                               log_prior_ratioNow +
                               log_proposal_ratioNow))
  
  if (runif(n = 1) < alpha_accept) {
    accepted = TRUE
    #log_print("Move ACCEPTED", console = FALSE)
    rho_updated = rho_proposed
  } else {
    accepted = FALSE
    #log_print("Move REJECTED", console = FALSE)
    rho_updated = rho_current
  }
  return(
    list(
      "rho_updated" = rho_updated,
      "accepted" = as.numeric(accepted),
      "choose_add" = choose_add,
      "candidate" = candidate,
      'alpha_accept' = alpha_accept,
      "three_ratios" = c('log_lik_ratio' = log_likelihood_ratioNow,
                         'log_prior_ratio' = log_prior_ratioNow,
                         'log_proposal_ratio' = log_proposal_ratioNow)
    )
  )
}



#' Adaptation step to update the weights vectors a and d
#'
#' The function takes as an input the current weights and updates them as a function of
#' the current iteration number t, the initial adaptation h
#' The function works in log scale
#'
#' @param logweights Vector of the logarithm of the current weights
#' @param alpha_target Scalar indicating the target acceptance probability (optimal range empirically observed around 0.10-0.15)
#' @param t Number of the current iteration
#' @param h Initial adaptation (must be >0)
#' @inheritParams update_partition
#'
#' @return Vector of updated logweights
#' @export
#'
#' @examples
#'
update_weight = function(weights, index, h, t_over_p, alpha_accept, alpha_target) {
    if (!h > 0)
        stop("Adaptation step h must be positive")
    # select the weight that has to be updated
    weight = weights[index]
    # update according to Benson
    weight = exp(log(weight) + h / t_over_p * (alpha_accept - alpha_target))
    # put it back
    weights[index] = weight
    # return the vector of weights
    return(weights)
}



#' Partition current data
#' Given y (vector of data) and rho (vector of the partition), the function splits the observations and
#' partitions them into the current groups
#' partitions them into the partition rho
#' @param y Vector of n ordered data
#' @inheritParams update_partition
#'
#' @return List where each element contains the corresponding group of y elements. If the dimensions of rho and y are not comparable, return an empty vector
#' @export
#'
#' @examples
partition_data <- function(y, rho) {
    if (sum(rho) != length(y))
        stop("The partition is not coherent with the data")
    
    dataPartition <- list()
    
    # number of groups
    M = length(rho)
    
    for (i in 1:M) {
        if (i == 1) {
            dataPartition[[i]] <- y[1:rho[i]]
            cumsum_rho = rho[i]
        } else {
            first_index = cumsum_rho + 1
            last_index = rho[i] + cumsum_rho
            dataPartition[[i]] <- y[first_index:last_index]
            cumsum_rho = cumsum_rho + rho[i]
        }
    }
    return(dataPartition)
}






#' Proposal Ratio
#'
#' @param choose_add Boolean to tell if we are performing an add move or not.
#'
#' @return The proposal ratio (not in log).
#' @export
#'
#' @examples
proposal_ratio = function(rho,
                          alpha_add,
                          weights_a,
                          weights_d,
                          choose_add) {
    # number of groups
    M = length(rho)
    
    n_elems = length(weights_a)
    
    # indexes of the changepoints
    cp_indexes <- get_group_indexes(rho)
    
    # exclude the last one because it's technically always 1 (a changepoint)
    cp_indexes <- cp_indexes[-length(cp_indexes)]
    
    # not all points can be selected for an add move
    # assign probability zero to those who cannot be
    weights_a_available = weights_a
    weights_a_available[cp_indexes] = 0
    weights_a_available_sum = sum(weights_a_available)
    
    # not all points can be selected for a delete move
    # assign probability zero to those who cannot be
    weights_d_available = weights_d
    weights_d_available[-cp_indexes] = 0
    weights_d_available_sum = sum(weights_d_available)
    
    if (choose_add) {
        draw_weights = weights_a_available
    } else {
        draw_weights = weights_d_available
    }
    
    # draw the candidate among the first 1:(p-1)
    candidate = sample(1:n_elems, 1, prob = draw_weights)
    
    if (choose_add && M == 1) {
        # case in which you choose to propose an add move
        # (with just 1 group) that may or may not be accepted
        ratio = (alpha_add / 1) * (weights_a_available_sum / weights_a[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    if (!choose_add && M == n_elems) {
        # case in which you choose to propose an delete move
        # (with every point being a group) that may or may not be accepted
        ratio = (alpha_add / 1) * (weights_d_available_sum / weights_d[candidate])
        return(list("ratio" = ratio, "candidate" = candidate))
    }
    
    # only the general cases remain
    if (choose_add) {
        ratio = (1 - alpha_add) / alpha_add *
            weights_a_available_sum / weights_a[candidate] *
            weights_d[candidate] / (weights_d[candidate] + weights_d_available_sum)
    } else {
        ratio = alpha_add / (1 - alpha_add) *
            weights_d_available_sum / weights_d[candidate] *
            weights_a[candidate] / (weights_a[candidate] + weights_a_available_sum)
    }
    
    return(list("ratio" = ratio, "candidate" = candidate))
}



#' Split Partition
#'
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the group that has changed.
#' @export
#'
#' @examples
split_partition <- function(candidate_index, rho) {

    # number of groups
    M = length(rho)
    new_rho = rep(NA, M + 1)
    
    group_indexes = get_group_indexes(rho)
    found = FALSE
    
    for (i in 1:M) {
        
        # update the partition in the general case
        # (either I have already split the group or not, just the index changes)
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i + 1] = rho[i]
        }
        
        if (!found && group_indexes[i] > candidate_index) {
            # just passed the element index - I am in the group to be split
            
            # index of the element minus the cumulative
            # number of elements in the previous groups only if i!=1
            new_rho[i] = candidate_index - (i != 1) * group_indexes[i - 1 * (i != 1)]
            
            # dimension of the original group minus the elements moved to new_rho[i]
            new_rho[i + 1] = rho[i] - new_rho[i]
            
            # save the index of the group that has changed
            j = i
            
            found = TRUE
        }
    }
    return(list("new_rho" = new_rho, "changed_group_index" = j))
}



#' Merge Partition
#'
#' @param candidate_index Index of the the point where to split the group (equivalent to adding a changepoint).
#' @inheritParams update_partition
#' @return A list whose first element is the updated partition and the second is the index of the group that has changed.
#' @export
#'
#' @examples
merge_partition <- function(candidate_index, rho) {

    # number of groups
    M = length(rho)
    new_rho = rep(NA, M - 1)
    
    group_indexes = get_group_indexes(rho)
    found = FALSE
    
    for (i in 1:(M - 1)) {
        
        # update the partition in the general case
        # either I have already merged the group or not, just the index changes
        if (!found) {
            new_rho[i] = rho[i]
        } else {
            new_rho[i] = rho[i + 1]
        }
        
        if (!found && group_indexes[i] == candidate_index) {
            # I am at the changepoint between the two groups to be merged
            
            # index of the element minus the cumulative
            # number of elements in the previous groups
            new_rho[i] = rho[i] + rho[i + 1]
            
            # save the index of the group that has changed
            j = i
            
            found = TRUE
        }
    }
    return(list("new_rho" = new_rho, "changed_group_index" = j))
}



#' Shuffle Partition
#'
#' @param G Adjacency matrix of the Graph.
#' @inheritParams update_partition
#' @return The updated partition after a shuffle move, if accepted.
#' @export
#'
#' @examples
shuffle_partition <- function(rho_current, G, sigma_prior, alpha, beta) {
    # vedi Corradin p.16 step (ii)
    
    # number of groups
    M = length(rho_current)
    
    # shuffling can be done only if the number of groups is at least 2
    if (M < 2) return(rho_current)
    
    # build new proposed rho
    rho_proposed = rho_current

    # going to shuffle group K with group K+1
    K <- sample(1:(M - 1), 1)
    
    if (rho_current[K] == 1 && rho_current[K+1] == 1){
        #log_print("SHUFFLE: cannot shuffle anything without reproposing the same partition", console = FALSE)
        return(rho_current)
    }
    
    # sample how many elements to keep in the K-th group
    sample_vector <- 1:(rho_current[K] + rho_current[K + 1] - 1)
    # avoid proposing the same partition (shuffling nothing)
    sample_vector <- sample_vector[-c(rho_current[K])]
    l <- sample(sample_vector, 1)

    # move the elements
    rho_proposed[K + 1] <- rho_current[K + 1] + rho_current[K] - l
    rho_proposed[K] <- l
    
    # compute log_prior_ratio
    log_prior_ratio = lpochhammer(1 - sigma_prior, rho_proposed[K])
                    + lpochhammer(1 - sigma_prior, rho_proposed[K + 1])
                    - lpochhammer(1 - sigma_prior, rho_current[K] - 1)
                    - lpochhammer(1 - sigma_prior, rho_current[K + 1] - 1)
                    + lfactorial(rho_current[K])
                    + lfactorial(rho_current[K + 1])
                    - lfactorial(rho_proposed[K])
                    - lfactorial(rho_proposed[K + 1])
    
    # to compute log_likelihood_ratio I need S

    S_current = get_S_from_G_rho(G, rho_current)
    S_star_current = get_S_star_from_S_and_rho(S_current, rho_current)
    
    S_proposed = get_S_from_G_rho(G, rho_proposed)
    S_star_proposed = get_S_star_from_S_and_rho(S_proposed, rho_proposed)
    
    # wrap the general fB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    fB = function(group1, group2, S, S_star) {
        return(fB_general(
            group1,
            group2,
            S,
            S_star,
            alpha = alpha,
            beta = beta,
            log = TRUE
        ))
    }
    
    # here there's not the ratio with the Beta coefficients
    # compute log_likelihood_ratio
    log_likelihood_ratio = 0
    
    for (l in 1:(K - 1)) {
        if (l > (K - 1)) break; # needed because R for-loops suck
        # first numerator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_current, S_star_current)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_current, S_star_current)
    }
    
    for (l in (K + 2):M) {
        if (l > M) break; # needed because R for-loops suck
        # second numerator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # second denominator term
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K, S_current, S_star_current)
        log_likelihood_ratio = log_likelihood_ratio + fB(l, K + 1, S_current, S_star_current)
    }
    
    # third numerator term
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K + 1, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K, S_proposed, S_star_proposed)
    log_likelihood_ratio = log_likelihood_ratio + fB(K + 1, K + 1, S_proposed, S_star_proposed)
    # third denominator term
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K + 1, S_current, S_star_current)
    log_likelihood_ratio = log_likelihood_ratio + fB(K, K, S_current, S_star_current)
    log_likelihood_ratio = log_likelihood_ratio + fB(K + 1, K + 1, S_current, S_star_current)

    # compute alpha_shuffle
    alpha_shuffle = min(1, exp(log_likelihood_ratio + log_prior_ratio))
    
    #log_print("SHUFFLE proposal", console = FALSE)
    #log_print("rho_current", console = FALSE)
    #log_print(rho_current, console = FALSE)
    #log_print("rho_proposed", console = FALSE)
    #log_print(rho_proposed, console = FALSE)

    if (runif(n = 1) < alpha_shuffle) {
        # accept the shuffle
        #log_print("Shuffle ACCEPTED", console = FALSE)
        return(rho_proposed)
    } else {
        # reject the shuffle
        #log_print("Shuffle REJECTED", console = FALSE)
        return(rho_current)
    }
}






# Build the entire S (sum of edges between clusters) from scratch
# idea: extracting all submatrices needed from G and sum all the elements
# (which are all ones). Do this only for the triangular part, then make it
# symmetric. You need to know just G and the partition rho
get_S_from_G_rho = function(G, rho) {
    
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups
    M = length(rho)
    
    # initialize S matrix
    S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
    
    # indexes of the right bounds of the partition
    bounds = cumsum(rho)
    
    # loop through the groups
    for (l in 1:M) {
        # extract the submatrix and sum all the elements
        for (m in 1:l) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            
            if(l == m){
                # the inside connections are now counted twice, correct for it
                S[l, m] = S[l, m] / 2
            } else {
                # otherwise, write to symmetric part of the matrix as well
                S[m, l] = S[l, m]
            }
        }
    }
    return(S)
}

# Build S (sum of edges between clusters) from knowledge of the previous S,
# the G matrix, the previous rho and the new proposed rho.
# Handles all the cases: add, delete, shuffle, i.e. "a new group is added",
# "a group is deleted", "same number of groups, but two groups exchange some
# elements".
get_S_from_G_rho_oldrho_oldS = function(G,rho,oldrho,oldS){
    
    if (!all(t(G) == G))
        stop("G is not symmetric")
    
    # number of groups in new and old rho
    M    = length(   rho)
    oldM = length(oldrho)
    
    # indexes of the right bounds of the partition
    bounds = get_group_indexes(rho)
    
    # groups that needs to be updated with the new rho information
    groups_to_be_refilled = {}
    
    if(M > oldM){ # case Add
        
        # initialize S matrix
        S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
        
        # find the group that has changed
        K = min(which(rho != c(oldrho,NA)))
        
        if(K > 1 && K < oldM){ # in the standard case perform all four
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
            # lower left block
            S[(K+1+1):M,1:(K-1)] = oldS[(K+1):oldM,1:(K-1)]
            
            # upper right block
            S[1:(K-1),(K+1+1):M] = oldS[1:(K-1),(K+1):oldM]
            
            # lower right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM]
            
        } else if(K == 1){ # bring to the new S only the lower right block
            
            # lower right block
            S[(K+1+1):M,(K+1+1):M] = oldS[(K+1):oldM,(K+1):oldM]
            
        } else if(K == oldM){ # bring to the new S only the upper left block
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
        }
        
        groups_to_be_refilled = c(K,K+1)
        
    } else if(M < oldM) { # case Delete
        
        # initialize S matrix
        S = matrix(numeric(M * M), nrow = M, byrow = TRUE)
        
        # find the group that has changed
        K = min(which(oldrho != c(rho,NA)))
        
        if(K > 1 && K+1 < oldM){ # in the standard case perform all four
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
            # lower left block
            S[(K+1):M,1:(K-1)] = oldS[(K+1+1):oldM,1:(K-1)]
            
            # upper right block
            S[1:(K-1),(K+1):M] = oldS[1:(K-1),(K+1+1):oldM]
            
            # lower right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM]
            
        } else if(K == 1){ # bring to the new S only the lower right block
            
            # lower right block
            S[(K+1):M,(K+1):M] = oldS[(K+1+1):oldM,(K+1+1):oldM]
            
        } else if(K+1 == oldM){ # bring to the new S only the upper left block
            
            # upper left block
            S[1:(K-1),1:(K-1)] = oldS[1:(K-1),1:(K-1)]
            
        }
        
        groups_to_be_refilled = c(K)
        
    } else { # case Shuffle
        
        S = oldS
        
        # find the group that has changed
        K = min(which(oldrho != rho))
        
        # set to zero the columns and rows of the shuffled groups
        S[K:(K+1),1:M] = 0 # row K to K+1
        S[1:M,K:(K+1)] = 0 # column K to K+1
        
        groups_to_be_refilled = c(K,K+1)
        
    }
    
    # loop through the groups
    for (l in groups_to_be_refilled) {
        # extract the submatrix and sum all the elements
        for (m in 1:M) {
            start_row = ifelse(m != 1, bounds[m - 1] + 1, 0)
            end_row = bounds[m]
            start_col = ifelse(l != 1, bounds[l - 1] + 1, 0)
            end_col = bounds[l]
            
            S[l, m] = sum(G[start_row:end_row, start_col:end_col])
            if(l == m){
                # the inside connections are now counted twice, correct for it
                S[l, m] = S[l, m] / 2
            } else {
                # otherwise, write to symmetric part of the matrix as well
                S[m, l] = S[l, m]
            }
        }
    }
    
    return(S)
}


#' Get Non-Edges from S and partition
#'
#' @param S MxM matrix, where M is number of groups, containing the sum of edges.
#' @inheritParams update_partition
#'
#' @return a positive scalar indicating the number of non-edges between two groups
#' computed as the possible number of edges between two groups
#' (depending on group cardinality) minus the effective number of edges
#' (depending on the edges actually present in the current Graph)
#' @export
#'
#' @examples
get_S_star_from_S_and_rho = function(S, rho){
    # number of groups
    M = length(rho)
    
    # initialize S_star matrix
    S_star = matrix(numeric(M * M), nrow = M, byrow = TRUE)
    
    # loop through the groups
    for (l in 1:M) {
        for (m in 1:l) {
            if (l == m){
                S_star[l,m] = rho[l] * (rho[l] - 1) / 2 - S[l,m]
            } else {
                S_star[l,m] = rho[l] * rho[m] - S[l,m]
                S_star[m,l] = S_star[l,m]
            }
        }
    }
    
    return(S_star)
}


# auxiliary function to evaluate the beta function for the likelihood ratio
fB_general = function(group1,
                      group2,
                      S,
                      S_star,
                      alpha,
                      beta,
                      log = TRUE) {
    if (log) {
        return(lbeta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    } else{
        return(beta(alpha + S[group1, group2], beta + S_star[group1, group2]))
    }
}

# auxiliary function to evaluate the beta function for the likelihood ratio
fB_zero = function(alpha,
                   beta,
                   log = TRUE) {
    if (log) {
        return(lbeta(alpha, beta))
    } else{
        return(beta(alpha, beta))
    }
}


#' Log Likelihood Ratio
#'
#' @inheritParams update_partition
#' @inheritParams shuffle_partition
#' @param rho_current Current partition.
#' @param rho_proposed Proposed partition.
#' @param alpha Parameter of the Beta of the Graph.
#' @param beta Parameter of the Beta of the Graph.
#'
#' @return The likelihood ratio in log.
#' @export
#'
#' @examples
log_likelihood_ratio = function(alpha_add,
                                weights_a,
                                weights_d,
                                G,
                                rho_current,
                                rho_proposed,
                                choose_add,
                                alpha,
                                beta,
                                changed_group_index) {
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = rho_current
        rho_current = rho_proposed
        rho_proposed = temp
    }
    
    # number of groups
    M = length(rho_current)
    
    # wrap the general fB into this version where I don't have to specify
    # alpha and beta (doing this for adaptiveness)
    fB = function(group1, group2, S, S_star) {
        return(fB_general(
            group1,
            group2,
            S,
            S_star,
            alpha = alpha,
            beta = beta,
            log = TRUE
        ))
    }
    
    S_current = get_S_from_G_rho(G, rho_current)
    S_star_current = get_S_star_from_S_and_rho(S_current, rho_current)
    
    S_proposed = get_S_from_G_rho(G, rho_proposed)
    S_star_proposed = get_S_star_from_S_and_rho(S_proposed, rho_proposed)
    
    K = changed_group_index
    #K = get_index_changed_group(rho_current, rho_proposed)
    
    # TODO
    # anziché tenere questo termine, metterlo in ognuno di quelli sotto.
    # Inefficiente, ma serve per consentire l'aggiornamento di alpha e beta
    log_ratio = -(M + 1) * fB_zero(alpha, beta)
    
    for (l in 1:(K - 1)) {
        if (l > (K - 1)) break; # needed because R for-loops suck
        # first numerator term
        log_ratio = log_ratio + fB(l, K, S_proposed, S_star_proposed)
        log_ratio = log_ratio + fB(l, K + 1, S_proposed, S_star_proposed)
        # first denominator term
        log_ratio = log_ratio - fB(l, K, S_current, S_star_current)
    }
    
    for (m in (K + 2):(M + 1)) {
        if (m > (M + 1)) break; # needed because R for-loops suck
        # second numerator term
        log_ratio = log_ratio + fB(K, m, S_proposed, S_star_proposed) +
                                fB(K + 1, m, S_proposed, S_star_proposed)
    }
    
    # third numerator term
    log_ratio = log_ratio + fB(K, K + 1, S_proposed, S_star_proposed) +
                            fB(K, K, S_proposed, S_star_proposed) +
                            fB(K + 1, K + 1, S_proposed, S_star_proposed)
    
    for (m in (K + 1):M) {
        if (m > M) break; # needed because R for-loops suck
        # second denominator term
        log_ratio = log_ratio - fB(K, m, S_current, S_star_current)
    }
    
    # third denominator term
    log_ratio = log_ratio - fB(K, K, S_current, S_star_current)
    
    # in the delete/merge case we have to invert everything
    if (!choose_add) {
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
    
}


#' Get index of the changed group between two partitions
#'
#' @param rho_current Current partition in the form of group cardinalities.
#' @param rho_proposed Proposed partition in the form of group cardinalities.
#'
#' @return Index of the group that has been affected by a change (split/merge/shuffle) from the current patition to the proposed one.
#' @export
#'
#' @examples
get_index_changed_group = function(rho_current, rho_proposed) {
    
    # indexes of the changepoints in the current partition
    cp_idxs_current = get_group_indexes(rho_current)
    
    # move from the rho representation to r representation
    # i.e. from c(2,3) to c(0,1,0,0,1)
    # both for the current and the proposed partition
     
    current_r = rho_to_r(rho_current)
    proposed_r = rho_to_r(rho_proposed)
    
    # now by making the difference I can extract the index
    # of the new changepoint (or the one deleted)
    # there's no need for an absolute value because in the delete move
    # they have been fictitiously swapped to avoid repeating code
    tau = which.max(abs(proposed_r - current_r))
    
    # take all the indexes of the cp smaller than tau
    temp = which(cp_idxs_current < tau)
    if(length(temp) == 0){
        temp = c(0)
    }
    
    # take the last element of this list to have the
    # index of the group affected
    K = temp[length(temp)] + 1
    
    return(K)
}


# TODO Add documentation
log_prior_ratio = function(theta_prior,
                          sigma_prior,
                          rho_current,
                          rho_proposed,
                          choose_add,
                          changed_group_index)
{
    # differentiate delete/merge case
    if (!choose_add) {
        # swap rhos 'cause we're lazy
        temp = rho_current
        rho_current = rho_proposed
        rho_proposed = temp
    }
    
    # number of groups in the current partition
    M = length(rho_current)
    
    K = changed_group_index
    #K = get_index_changed_group(rho_current,rho_proposed)
    
    # compute the prior ratio
    log_ratio = - log(M) + log(theta_prior + M * sigma_prior)
                + lpochhammer(1 - sigma_prior, rho_proposed[K] - 1)
                + lpochhammer(1 - sigma_prior, rho_proposed[K + 1] - 1)
                - lpochhammer(1 - sigma_prior, rho_proposed[K] + rho_proposed[K + 1] - 1)
                + lfactorial(rho_proposed[K] + rho_proposed[K + 1])
                - lfactorial(rho_proposed[K])
                - lfactorial(rho_proposed[K + 1])
    
    
    # in the delete/merge case we have to invert everything
    if (!choose_add) {
        log_ratio = -log_ratio
    }
    
    return(log_ratio)
}




#' Compute weights of the full conditional of theta
#' 
#' For further details see Proposition 1 Martinez and Mena (2014).
#'
#' @param p number of nodes
#' @param sigma_prior other parameter used to compute the prior ratio
#' @param k number of groups
#' @param j index of the iteration for which we are computing the weight
#' @param f value drawn from Exp(theta+1)
#' @param z value drawn from Be(theta+2,n)
#' @param c prior first parameter of the shifted gamma
#' @param d prior second parameter of the shifted gamma
#'
#' @return
#' @export
#'
#' @examples
compute_weights_theta <- function(c, d, p, sigma_prior, k, j, f, z) {
    abs_stir <- abs_stirling_number_1st(k, j)
    num <- (
        (p - sigma_prior) * (p + 1 - sigma_prior) * abs_stirling_number_1st(k - 1, j) +
            (2 * p + 1 - 2 * sigma_prior) * sigma_prior * abs_stirling_number_1st(k - 1, j - 1) +
            (sigma_prior ^ 2) * abs_stirling_number_1st(k - 1, j - 2)
    ) * gamma(c + j)

    denom <- (sigma_prior * (d + f - log(z)))^j
    return(num / denom)
}



#' Full-conditional for theta
#' 
#' For further details see Proposition 1 Martinez and Mena (2014).
#' 
#' TODO COMPLETE DOCUMENTATION
#'
#' @param c First parameter of the shifted gamma prior
#' @param d Second parameter of the shifted gamma prior
#' @param candidate proposed value for theta
#' @param k number of groups
#'
#' @return scalar value for theta at the current iteration
#' @export
#'
#' @examples
full_conditional_theta <- function(prior_c, prior_d, candidate, k, p, sigma_prior){
    weights_gamma <- rep(0,k+2)
    z = rbeta(1,candidate + 2, p)
    f = rexp (1,candidate + 1)
    
    for (j in 0:(k+1)){
        weight_j = compute_weights_theta(prior_c, prior_d, p, sigma_prior, k, j, f, z)
        weights_gamma[j] = weight_j
    }

    # normalizing the weights
    if(sum(weights_gamma) != 0) {
        weights_gamma = weights_gamma / sum(weights_gamma)
    }
    
    component = min(which(cumsum(weights_gamma) > runif(n = 1)))
    
    theta = shifted_gamma(prior_c + (component - 1), prior_d + f - log(z), -sigma_prior)
    return(theta)
}



#' Full-conditional for sigma
#' 
#' For further details see Section 4 Martinez and Mena (2014).
#'
#' @param sigma value of sigma to be updated
#' @param theta value of theta
#' @param rho Current partition
#' @param a First parameter of the distribution
#' @param b Second parameter of the distribution
#' @param c Third parameter of the distribution
#' @param d Fourth parameter of the distribution
#'
#' @return the value for sigma for the current iteration
#' @export
#'
#' @examples
full_conditional_sigma <- function(sigma, theta, rho, a, b, c, d){
    
    # number of groups
    M <- length(rho)

    # first product term
    log_prod_1 <- 0

    for (i in 1:(M - 1)) {
        if (i > (M - 1)) break; # needed because R for-loops suck
        log_prod_1 <- log_prod_1 + log(theta + i * sigma)
    }

    # second product term
    log_prod_2 <- 0
    for (i in 1:M) {
        log_prod_2 <- log_prod_2 + lpochhammer((1 - sigma), (rho[i] - 1))
    }

    # final output
    output <- (a - 1) * log(sigma) + (b - 1) * log(1 - sigma) +
        (c - 1) * log(theta + sigma) + (-d * sigma) +
        log_prod_1 + log_prod_2

    return(output)
}


#' Set options for the Gibbs sampler
#'
#' This function sets various options required for running the Gibbs sampler.
#' It allows customization of parameters, such as priors, initial values, adaptation steps, update flags, etc.
#' 
#' @param sigma_prior_0 Initial parameter of the prior (sigma) from Martinez and Mena (2014).
#' @param sigma_prior_parameters List of parameters for updating sigma_prior.
#' @param theta_prior_0 Initial parameter of the prior (theta) from Martinez and Mena (2014).
#' @param theta_prior_parameters List of parameters for updating theta_prior.
#' @param rho0 Initial partition.
#' @param weights_a0 Weights for choosing the candidate node in an add/split move.
#' @param weights_d0 Weights for choosing the candidate node in a delete/merge move.
#' @param total_weights0 Total weights (added parameter).
#' @param total_K0 Total K (added parameter).
#' @param total_graphs0 Total graphs (added parameter).
#' @param graph Initial graph.
#' @param alpha_target Target acceptance rate of the split and merge Metropolis-Hastings.
#' @param beta_mu Expected value for the Beta prior of the graph.
#' @param beta_sig2 Variance for the Beta prior of the graph.
#' @param d Parameter of the G-Wishart.
#' @param alpha_add Probability of choosing an add/split move over a delete/merge move.
#' @param adaptation_step Adaptation step for tweaking how much the weights are updated each time.
#' @param update_sigma_prior Boolean for choosing whether to update sigma_prior or not.
#' @param update_theta_prior Boolean for choosing whether to update theta_prior or not.
#' @param update_weights Boolean for choosing whether to update weights or not.
#' @param update_partition Boolean for choosing whether to update the partition or not.
#' @param update_graph Boolean for choosing whether to update the graph or not.
#' @param perform_shuffle Boolean for choosing whether to perform shuffle or not.
#' 
#' @return A list with all options correctly set for working with the Gibbs sampler.
#' @export
#'
#' @examples
set_options = function(sigma_prior_0,
                       sigma_prior_parameters,
                       theta_prior_0,
                       theta_prior_parameters,
                       rho0,
                       weights_a0,
                       weights_d0,
                       total_weights0, 
                       total_K0, 
                       total_graphs0, 
                       graph,   
                       alpha_target,
                       beta_mu,
                       beta_sig2,
                       d=3,
                       alpha_add=0.5,
                       adaptation_step,
                       update_sigma_prior=TRUE,
                       update_theta_prior=TRUE,
                       update_weights=TRUE,
                       update_partition=TRUE,
                       update_graph=TRUE,
                       perform_shuffle=TRUE
) {
    
    options = list(
        "sigma_prior_0"          = sigma_prior_0,
        "sigma_prior_parameters" = sigma_prior_parameters,
        "theta_prior_0"          = theta_prior_0,
        "theta_prior_parameters" = theta_prior_parameters,
        "rho0"                   = rho0,
        "weights_a0"             = weights_a0,
        "weights_d0"             = weights_d0,
        "total_weights0"         = total_weights0,
        "total_K0"               = total_K0,
        "total_graphs0"          = total_graphs0,
        "graph"                  = graph,
        "alpha_target"           = alpha_target,
        "beta_mu"                = beta_mu,
        "beta_sig2"              = beta_sig2,
        "d"                      = d,
        "alpha_add"              = alpha_add,
        "adaptation_step"        = adaptation_step,
        "update_sigma_prior"     = update_sigma_prior,
        "update_theta_prior"     = update_theta_prior,
        "update_weights"         = update_weights,
        "update_partition"       = update_partition,
        "update_graph"           = update_graph,
        "perform_shuffle"        = perform_shuffle
    )
    return(options)
}


#' Estimate alpha and beta from a Beta function given mean and variance
#' See https://stats.stackexchange.com/a/12239 for details.
#'
#' @param mu Mean of the Beta.
#' @param var Variance of the Beta.
#'
#' @return List with alpha and beta of the corresponding Beta distribution.
#' @export
#'
#' @examples
estimate_Beta_params <- function(mu, var) {
    if(!(var > 0 && var < mu*(1-mu)))
        stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")
    
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(list(alpha = alpha, beta = beta))
}


#' Gibbs sampler
#'
#' This function implements the Gibbs sampler algorithm for Bayesian inference.
#' It iteratively samples from the posterior distribution of parameters using conditional distributions.
#' 
#' @param data An n x p matrix of data.
#' @param niter Desired number of effective iterations.
#' @param nburn Number of iterations to be burned.
#' @param thin Keep all multiples of thin.
#' @param options All the parameters necessary to run the Gibbs sampler.
#' @param seed Seed for reproducibility.
#' @param print Print the progress bar. Default to TRUE.
#'
#' @return A list containing sampled values of parameters and other information.
#' @export
#'
#' @examples
# Gibbs_sampler = function(data,
#                          niter,
#                          nburn,
#                          thin,
#                          options,
#                          seed=1234,
#                          print=TRUE,
#                          s)
# {
#     n = nrow(data) # number of observations
#     p = ncol(data) # number of nodes
#     n_total_iter = nburn + niter * thin # total iterations to be made
#     
#     # dynamic parameters
#     sigma_prior            = options$sigma_prior_0
#     theta_prior            = options$theta_prior_0
#     rho                    = options$rho0
#     weights_a              = options$weights_a0
#     weights_d              = options$weights_d0
#     total_weights          = options$total_weights0     # added
#     total_K                = options$total_K0           # added
#     total_graphs           = options$total_graphs0      # added
#     graph                  = options$graph              # graph start or bd_graph_start
#     last_S                 = NULL
#     adaptation_step        = options$adaptation_step
#     sigma_prior_parameters = options$sigma_prior_parameters
#     theta_prior_parameters = options$theta_prior_parameters
#     
#     # constant parameters
#     alpha_add    = options$alpha_add
#     alpha_target = options$alpha_target
#     
#     # parameter for the Wishart
#     d = options$d
#     
#     # parameters for the Beta
#     beta_mu   = options$beta_mu
#     beta_sig2 = options$beta_sig2
#     
#     
#     # checks
#     if(sum(rho) != p)
#         stop("The partition rho must sum to the number of variables p")
#     if(d < 3)
#         stop("The Wishart's d parameter must be greater or equal than 3")
#     if(!(beta_mu > 0 && beta_mu < 1))
#         stop("The mean of the Beta must be between 0 and 1")
#     if(!(beta_sig2 > 0 && beta_sig2 < beta_mu*(1-beta_mu)))
#         stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")
#     if(length(weights_a) != (p-1) || length(weights_d) != (p-1))
#         stop("The number of elements in the weights vectors must be equal to p-1")
#     if(!(adaptation_step > 0))
#         stop("The adapation step h must be positive")
#     if(!(alpha_add > 0 && alpha_add < 1))
#         stop("The probability of choosing an add move alpha_add must be between 0 and 1")
#     if(!(alpha_target > 0 && alpha_target < 1))
#         stop("The target acceptance rate of the Metropolis-Hastings alpha_target must be between 0 and 1")
#     
#     # t_over_p = n_total_iter / p
#     t_over_p = s / p #Luciano: this for adaptivity
#     
#     beta_params = estimate_Beta_params(beta_mu, beta_sig2)
#     
#     # define structure to save sampled values
#     save_res = list(
#         G = vector("list", length = niter),
#         K = vector("list", length = niter),
#         rho = vector("list", length = niter),
#         accepted = vector("numeric", length = niter),
#         S = vector("list", length = niter),
#         sigma = vector("numeric", length = niter),
#         theta = vector("numeric", length = niter),
#         weights_a = vector("list", length = niter),
#         weights_d = vector("list", length = niter),
#         bdgraph_start = NULL,
#         total_graphs = vector("list", length = niter),
#         total_K = vector("list", length = niter),
#         total_weights = 0
#     )
#     
#     # initialize iteration counter
#     it_saved = 0
#     
#     # initialize progress bar
#     if(print){
#         pb = txtProgressBar(min=1, max=n_total_iter, initial=1, style=3)
#     }
#     
#     # save start time for measuring execution time
#     start_time = Sys.time()
#     
#     # start the simulation
#     for(iter in 1:n_total_iter){
#         
#         if(is.null(graph)){ 
#             # we run a single iteration of BDgraph with iter = 1 and burnin = 0
#             # and g.start = "empty"
#             output = bdgraph(
#                 data,
#                 rho,
#                 n,
#                 method = "ggm",
#                 algorithm = "bdmcmc",
#                 iter = 1,
#                 burnin = 0,
#                 not.cont = NULL,
#                 g.prior = 0.5,
#                 df.prior = d,
#                 CCG_D = NULL,
#                 g.start = "empty",
#                 jump = NULL,
#                 save = TRUE,
#                 print = 1000,
#                 cores = NULL,
#                 threshold = 1e-8
#             )
#         }
#         
#         else {
#             # we run a single iteration of BDgraph with iter = 1 and burnin = 0
#             # and g.start = graph
#             output = bdgraph(
#                 data,
#                 rho,
#                 n,
#                 method = "ggm",
#                 algorithm = "bdmcmc",
#                 iter = 1,
#                 burnin = 0,
#                 not.cont = NULL,
#                 g.prior = 0.5,
#                 df.prior = d,
#                 CCG_D = NULL,
#                 g.start = graph,
#                 jump = NULL,
#                 save = TRUE,
#                 print = 1000,
#                 cores = NULL,
#                 threshold = 1e-8
#             )
#         }
#         
#         # update graph
#         if (options$update_graph){
#             
#             # extract precision matrix K
#             last_K = output$last_K
#             
#             if(niter > nburn){ # only if niter > nburn right? TODO what about burnin? Siamo 'sgor?
#                 # update total_weights
#                 total_weights = total_weights + output$all_weights
#                 # update total_graphs taking into consideration the weights
#                 total_graphs = total_graphs + output$last_graph * output$all_weights
#                 total_K = total_K + output$last_K * output$all_weights
#             }
#             
#         }else{
#           # extract precision matrix K even if not updating graph
#           last_K = output$last_K  
#         }
#         
#         # extract adjacency matrix G and update for the next iteration
#         graph = output$last_graph
#         # e mettere in update partition non graph ma total_graphs/weights
#         # se update graph oppure total graph se non update weights??
#         
#         if(options$update_partition){
#             list_output_update_partition = update_partition(rho,
#                                                             alpha_add,
#                                                             weights_a,
#                                                             weights_d,
#                                                             theta_prior,
#                                                             sigma_prior,
#                                                             graph,
#                                                             beta_params)
#             
#             rho = list_output_update_partition$rho_updated
#             
#             # it makes sense to perform the adaptive step only if we're updating the partition
#             # update the single weight at the point only if the move has been accepted
#             
#             
#             if(options$update_weights && list_output_update_partition$accepted){
#                 if(list_output_update_partition$choose_add){
#                     weights_a = update_weight(
#                         weights_a,
#                         list_output_update_partition$candidate,
#                         adaptation_step,
#                         t_over_p,
#                         alpha_add,
#                         alpha_target
#                     )
#                 } else {
#                     weights_d = update_weight(
#                         weights_d,
#                         list_output_update_partition$candidate,
#                         adaptation_step,
#                         t_over_p,
#                         1 - alpha_add,
#                         alpha_target
#                     )
#                     
#                 }
#             }
#         }
#         
#         
#         if(options$perform_shuffle){
#             rho = shuffle_partition(rho, graph, sigma_prior, beta_params$alpha, beta_params$beta)
#         }
#         
#         if(options$update_sigma_prior){
#             candidate <- runif(1,max(0,-theta_prior),1)
#             alpha_MH <- full_conditional_sigma(candidate,
#                                                theta_prior,
#                                                rho,
#                                                sigma_prior_parameters$a,
#                                                sigma_prior_parameters$b,
#                                                sigma_prior_parameters$c,
#                                                sigma_prior_parameters$d) -
#                 full_conditional_sigma(sigma_prior,
#                                        theta_prior,
#                                        rho,
#                                        sigma_prior_parameters$a,
#                                        sigma_prior_parameters$b,
#                                        sigma_prior_parameters$c,
#                                        sigma_prior_parameters$d)
#             
#             if(log(runif(1)) <= min(alpha_MH,log(1))){
#                 sigma_prior = candidate
#             }else{
#                 sigma_prior = sigma_prior
#             }
#         }
#         
#         if(options$update_theta_prior) {
#             theta_prior = full_conditional_theta(
#                 theta_prior_parameters$c,
#                 theta_prior_parameters$d,
#                 theta_prior,
#                 length(rho),
#                 p,
#                 sigma_prior
#             )
#         }
#         
#         if(options$update_graph){
#             last_S = get_S_from_G_rho(graph,rho)
#         }
#         
#         # save results only on thin iterations
#         # (i.e. only save multiples of thin)
#         if(iter > nburn && (iter - nburn) %% thin == 0) {
#             it_saved = it_saved + 1
#             if(options$update_graph){
#                 # cumulative precision matrix K and probability of inclusion links
#                 save_res$K[[it_saved]] = total_K / total_weights
#                 save_res$G[[it_saved]] = total_graphs / total_weights
#             }
#             else{
#                 save_res$K[[it_saved]] = total_K
#                 save_res$G[[it_saved]] = total_graphs
#             }
#             
#             save_res$total_weights <- total_weights
#             save_res$total_K[[it_saved]] <- total_K
#             save_res$total_graphs[[it_saved]] <- total_graphs
#             
#             # save bdgraph object
#             # recall graph = output$last_graph is adjacency
#             save_res$bdgraph_start = output$last_graph
#             
#             if(options$update_partition){
#                 save_res$accepted[[it_saved]] = list_output_update_partition$accepted
#             }
#             
#             save_res$rho[[it_saved]] = rho
#             save_res$weights_d[[it_saved]] <- weights_d
#             save_res$weights_a[[it_saved]] <- weights_a
#             
#             save_res$sigma[[it_saved]] = sigma_prior
#             save_res$theta[[it_saved]] = theta_prior
#             save_res$S[[it_saved]] = last_S
#         }
#         
#         if(print){
#             setTxtProgressBar(pb, iter)
#         }
#         #log_print("last_G:", console = FALSE)
#         #log_print(last_G, console = FALSE)
#         #log_print("last_S:", console = FALSE)
#         #log_print(last_S, console = FALSE)
#         #log_print("---------------------------------------------------------------------", console = FALSE)
#     }
#     
#     save_res$execution_time = Sys.time() - start_time
#     
#     #output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = "empty",
#     #                   all_graphs = niter-nburn, all_weights = all_weights, last_graph = last_G,
#     #                   last_K = last_K, last_Theta = last_Theta )
#     #
#     
#     if(print){close(pb)}
#     return(save_res)
# }

Gibbs_sampler = function(data,
                         niter,
                         nburn,
                         thin,
                         options,
                         seed=1234,
                         print=TRUE,
                         s,
                         algorithm)
{
  n = nrow(data) # number of observations
  p = ncol(data) # number of nodes
  n_total_iter = nburn + niter * thin # total iterations to be made
  
  # dynamic parameters
  sigma_prior            = options$sigma_prior_0
  theta_prior            = options$theta_prior_0
  rho                    = options$rho0
  weights_a              = options$weights_a0
  weights_d              = options$weights_d0
  total_weights          = options$total_weights0     # added
  total_K                = options$total_K0           # added
  total_graphs           = options$total_graphs0      # added
  graph                  = options$graph              # graph start or bd_graph_start
  last_S                 = NULL
  adaptation_step        = options$adaptation_step
  sigma_prior_parameters = options$sigma_prior_parameters
  theta_prior_parameters = options$theta_prior_parameters
  
  # constant parameters
  alpha_add    = options$alpha_add
  alpha_target = options$alpha_target
  
  # parameter for the Wishart
  d = options$d
  
  # parameters for the Beta
  beta_mu   = options$beta_mu
  beta_sig2 = options$beta_sig2
  
  # checks
  if(sum(rho) != p)
    stop("The partition rho must sum to the number of variables p")
  if(d < 3)
    stop("The Wishart's d parameter must be greater or equal than 3")
  if(!(beta_mu > 0 && beta_mu < 1))
    stop("The mean of the Beta must be between 0 and 1")
  if(!(beta_sig2 > 0 && beta_sig2 < beta_mu*(1-beta_mu)))
    stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")
  if(length(weights_a) != (p-1) || length(weights_d) != (p-1))
    stop("The number of elements in the weights vectors must be equal to p-1")
  if(!(adaptation_step > 0))
    stop("The adapation step h must be positive")
  if(!(alpha_add > 0 && alpha_add < 1))
    stop("The probability of choosing an add move alpha_add must be between 0 and 1")
  if(!(alpha_target > 0 && alpha_target < 1))
    stop("The target acceptance rate of the Metropolis-Hastings alpha_target must be between 0 and 1")
  
  # t_over_p = n_total_iter / p
  t_over_p = s / p #Luciano: this for adaptivity
  
  beta_params = estimate_Beta_params(beta_mu, beta_sig2)
  
  # define structure to save sampled values
  save_res = list(
    G = vector("list", length = niter),
    K = vector("list", length = niter),
    rho = vector("list", length = niter),
    accepted = vector("numeric", length = niter),
    S = vector("list", length = niter),
    sigma = vector("numeric", length = niter),
    theta = vector("numeric", length = niter),
    weights_a = vector("list", length = niter),
    weights_d = vector("list", length = niter),
    bdgraph_start = NULL,
    total_graphs = vector("list", length = niter),
    total_K = vector("list", length = niter),
    total_weights = 0
  )
  
  # initialize iteration counter
  it_saved = 0
  
  # initialize progress bar
  if(print){
    pb = txtProgressBar(min=1, max=n_total_iter, initial=1, style=3)
  }
  
  # save start time for measuring execution time
  start_time = Sys.time()
  
  # start the simulation
  for(iter in 1:n_total_iter){
    
    if(is.null(graph)){ 
      # we run a single iteration of BDgraph with iter = 1 and burnin = 0
      # and g.start = "empty"
      output = post_graph_sampling(
        data,
        rho,
        n,
        method = "ggm",
        algorithm = algorithm,
        iter = 1,
        burnin = 0,
        not.cont = NULL,
        g.prior = 0.5,
        df.prior = d,
        CCG_D = NULL,
        g.start = "empty",
        jump = NULL,
        save = TRUE,
        print = 1,
        cores = NULL,
        threshold = 1e-8,
        beta_params = beta_params
      )
    } else {
      # we run a single iteration of BDgraph with iter = 1 and burnin = 0
      # and g.start = graph
      output = post_graph_sampling(
        data,
        rho,
        n,
        method = "ggm",
        algorithm = algorithm,
        iter = 1,
        burnin = 0,
        not.cont = NULL,
        g.prior = 0.5,
        df.prior = d,
        CCG_D = NULL,
        g.start = graph,
        jump = NULL,
        save = TRUE,
        print = 1,
        cores = NULL,
        threshold = 1e-8,
        beta_params = beta_params
      )
    }
    
    # update graph
    if (options$update_graph){
      
      # extract precision matrix K
      last_K = output$last_K
      
      if(algorithm == 'bdmcmc'){ 
        # update total_weights
        total_weights = total_weights + output$all_weights
        # update total_graphs taking into consideration the weights
        total_graphs = total_graphs + output$last_graph * output$all_weights
        total_K = total_K + output$last_K * output$all_weights
      }else{
        # in this case we are using rjmcmc hence every graph and precision have weight 1 
        # if doing 1 iteration
        total_weights = total_weights + output$all_weights 
        total_graphs = total_graphs + output$last_graph
        total_K = total_K + output$last_K
      }
    }else{
      # extract precision matrix K even if not updating graph
      last_K = output$last_K  
    }
    
    # extract adjacency matrix G and update for the next iteration
    graph = output$last_graph
    
    if(options$update_partition){
      list_output_update_partition = update_partition(rho,
                                                      alpha_add,
                                                      weights_a,
                                                      weights_d,
                                                      theta_prior,
                                                      sigma_prior,
                                                      graph,
                                                      beta_params)
      
      rho = list_output_update_partition$rho_updated
      alpha_accept = list_output_update_partition$alpha_accept
      # it makes sense to perform the adaptive step only if we're updating the partition
      # update the single weight at the point only if the move has been accepted
      
      if(options$update_weights && list_output_update_partition$accepted){
        if(list_output_update_partition$choose_add){
          weights_a = update_weight(
            weights_a,
            list_output_update_partition$candidate,
            adaptation_step,
            t_over_p,
            alpha_accept,
            alpha_target
          )
        } else {
          weights_d = update_weight(
            weights_d,
            list_output_update_partition$candidate,
            adaptation_step,
            t_over_p,
            alpha_accept,
            alpha_target
          )
          
        }
      }
    }
    
    if(options$perform_shuffle){
      rho = shuffle_partition(rho, graph, sigma_prior, beta_params$alpha, beta_params$beta)
    }
    
    if(options$update_sigma_prior){
      candidate <- runif(1,max(0,-theta_prior),1)
      alpha_MH <- full_conditional_sigma(candidate,
                                         theta_prior,
                                         rho,
                                         sigma_prior_parameters$a,
                                         sigma_prior_parameters$b,
                                         sigma_prior_parameters$c,
                                         sigma_prior_parameters$d) -
        full_conditional_sigma(sigma_prior,
                               theta_prior,
                               rho,
                               sigma_prior_parameters$a,
                               sigma_prior_parameters$b,
                               sigma_prior_parameters$c,
                               sigma_prior_parameters$d)
      
      if(log(runif(1)) <= min(alpha_MH,log(1))){
        sigma_prior = candidate
      }else{
        sigma_prior = sigma_prior
      }
    }
    
    if(options$update_theta_prior) {
      theta_prior = full_conditional_theta(
        theta_prior_parameters$c,
        theta_prior_parameters$d,
        theta_prior,
        length(rho),
        p,
        sigma_prior
      )
    }
    
    if(options$update_graph){
      last_S = get_S_from_G_rho(graph,rho)
    }
    
    # save results only on thin iterations
    # (i.e. only save multiples of thin)
    if(iter > nburn && (iter - nburn) %% thin == 0) {
      it_saved = it_saved + 1
      if(options$update_graph){
        # cumulative precision matrix K and probability of inclusion links
        save_res$K[[it_saved]] = total_K / total_weights
        save_res$G[[it_saved]] = total_graphs / total_weights
      }
      else{
        save_res$K[[it_saved]] = total_K
        save_res$G[[it_saved]] = total_graphs
      }
      
      save_res$total_weights <- total_weights
      save_res$total_K[[it_saved]] <- total_K
      save_res$total_graphs[[it_saved]] <- total_graphs
      
      # save bdgraph object
      # recall graph = output$last_graph is adjacency
      save_res$bdgraph_start = output$last_graph
      
      if(options$update_partition){
        save_res$accepted[[it_saved]] = list_output_update_partition$accepted
      }
      
      save_res$rho[[it_saved]] = rho
      save_res$weights_d[[it_saved]] <- weights_d
      save_res$weights_a[[it_saved]] <- weights_a
      
      save_res$sigma[[it_saved]] = sigma_prior
      save_res$theta[[it_saved]] = theta_prior
      save_res$S[[it_saved]] = last_S
    }
    
    if(print){
      setTxtProgressBar(pb, iter)
    }
    #log_print("last_G:", console = FALSE)
    #log_print(last_G, console = FALSE)
    #log_print("last_S:", console = FALSE)
    #log_print(last_S, console = FALSE)
    #log_print("---------------------------------------------------------------------", console = FALSE)
  }
  
  save_res$execution_time = Sys.time() - start_time
  
  #output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = "empty",
  #                   all_graphs = niter-nburn, all_weights = all_weights, last_graph = last_G,
  #                   last_K = last_K, last_Theta = last_Theta )
  #
  
  if(print){close(pb)}
  return(save_res)
}



#' Set initialization values for the Gibbs sampler
#'
#' This function sets the initialization values for the Gibbs sampler.
#' 
#' @param Beta Initial value for Beta.
#' @param mu Initial value for mu.
#' @param tau_eps Initial value for tau_eps.
#' @param K Initial value for K.
#' @param G Initial value for G.
#' @param z Initial value for z.
#' @param rho Initial value for rho.
#' @param a_sigma Hyperparameter of sigma prior.
#' @param b_sigma Hyperparameter of sigma prior.
#' @param c_sigma Hyperparameter of sigma prior.
#' @param d_sigma Hyperparameter of sigma prior.
#' @param c_theta Hyperparameter of theta prior.
#' @param d_theta Hyperparameter of theta prior.
#' @param sigma Initial value for sigma.
#' @param theta Initial value for theta.
#' @param weights_a0 Initial value for weights_a0.
#' @param weights_d0 Initial value for weights_d0.
#' @param total_weights Initial value for total_weights.
#' @param total_K Initial value for total_K.
#' @param total_graphs Initial value for total_graphs.
#' @param graph_start Initial value for graph_start.
#' @param graph_density Initial value for graph_density.
#' @param beta_sig2 Initial value for beta_sig2.
#' @param d Initial value for d.
#'
#' @return A list containing the initialization values for the Gibbs sampler.
#'
#' @examples
#'
set_initialization = function(
        Beta,
        mu,
        tau_eps,
        K,
        G,
        z,
        rho,
        a_sigma,
        b_sigma,
        c_sigma,
        d_sigma,
        c_theta,
        d_theta,
        sigma,
        theta,
        weights_a0,
        weights_d0,
        total_weights,
        total_K,
        total_graphs,
        graph_start,
        graph_density,
        beta_sig2,
        d
    ){
    
    initialization_values = list(
        'Beta'            = Beta,
        'mu'              = mu,
        'tau_eps'         = tau_eps,
        'K'               = K,
        'G'               = G,
        'z'               = z,
        'rho'             = rho,
        'a_sigma'         = a_sigma,
        'b_sigma'         = b_sigma,
        'c_sigma'         = c_sigma,
        'd_sigma'         = d_sigma,
        'c_theta'         = c_theta,
        'd_theta'         = d_theta,
        'sigma'           = sigma,
        'theta'           = theta,
        'weights_a'       = weights_a0,
        'weights_d'       = weights_d0,
        'total_weights'   = total_weights,
        'total_K'         = total_K,
        'total_graphs'    = total_graphs,
        'graph_start'     = graph_start,
        'graph_density'   = graph_density,
        'beta_sig2'       = beta_sig2,
        'd'               = d
    )
    
    return(initialization_values)
}



#' Set update parameters for Gibbs sampler
#'
#' This function sets the update parameters for the Gibbs sampler.
#' 
#' @param tbase_base Base for tbase.
#' @param tbase_data Data for tbase.
#' @param Sdata Data for S.
#' @param a_tau_eps Value for a_tau_eps.
#' @param b_tau_eps Value for b_tau_eps.
#' @param sigma_mu Value for sigma_mu.
#' @param r Value for r.
#' @param Update_Beta Boolean for updating Beta.
#' @param Update_Mu Boolean for updating Mu.
#' @param Update_Tau Boolean for updating Tau.
#'
#' @return A list containing the update parameters.
#'
#' @examples
#'
set_UpdateParamsGSL = function(
    tbase_base,
    tbase_data,
    Sdata,
    a_tau_eps,
    b_tau_eps,
    sigma_mu,
    r,
    Update_Beta = TRUE,
    Update_Mu = TRUE,
    Update_Tau = TRUE
){
    return(
        list(
        'tbase_base'   = tbase_base,
        'tbase_data'   = tbase_data,
        'Sdata'        = Sdata,
        'a_tau_eps'    = a_tau_eps,
        'b_tau_eps'    = b_tau_eps,
        'sigma_mu'     = sigma_mu,
        'r'            = r,
        'Update_Beta'  = Update_Beta,
        'Update_Mu'    = Update_Mu,
        'Update_Tau'   = Update_Tau
        )
    )
}


#' Update the Gibbs sampler
#'
#' This function updates the Gibbs sampler for a specified number of iterations.
#' 
#' @param set_UpdateParamsGSL_list List of update parameters.
#' @param niter Number of iterations.
#' @param initialization_values Initial values for the Gibbs sampler.
#' @param alpha_target Target acceptance rate for Metropolis-Hastings.
#' @param alpha_add Probability of choosing an add move.
#' @param adaptation_step Step size for adaptation.
#' @param seed Seed for reproducibility.
#' @param update_sigma_prior Boolean for updating sigma prior.
#' @param update_theta_prior Boolean for updating theta prior.
#' @param update_weights Boolean for updating weights.
#' @param update_partition Boolean for updating partition.
#' @param update_graph Boolean for updating graph.
#' @param perform_shuffle Boolean for shuffling.
#'
#' @return A list containing the updated results of the Gibbs sampler.
#'
#' @examples
#'
# Gibbs_sampler_update = function(
#         set_UpdateParamsGSL_list,
#         niter, 
#         initialization_values, 
#         alpha_target,
#         alpha_add,
#         adaptation_step,
#         seed,
#         update_sigma_prior,
#         update_theta_prior, 
#         update_weights,
#         update_partition,
#         update_graph,
#         perform_shuffle
# ){
#     set.seed(seed)
#     # Create a list for chains to save the values of each iteration
#     chains <- list(
#         Beta = vector("list", length = niter),
#         mu = vector("list", length = niter),
#         tau_eps = vector("list", length = niter),
#         K = vector("list", length = niter),
#         G = vector("list", length = niter),
#         z = vector("list", length = niter),
#         rho = vector("list", length = niter),
#         time = vector("list", length = niter), 
#         sigma = vector("list", length = niter),
#         theta = vector("list", length = niter),
#         weights_a = vector("list", length = niter),
#         weights_d = vector("list", length = niter),
#         # da togliere eventualmente
#         bd_graph_start = vector("list", length = niter), # to check graph G for likelihood of partition
#         total_weights = vector("list", length = niter) # to check graph G for likelihood of partition
#     ) 
#     
#     # Initialization of the chains
#     chains$Beta[[1]] <- initialization_values$Beta        
#     chains$mu[[1]] <- initialization_values$mu
#     chains$tau_eps[[1]] <- initialization_values$tau_eps 
#     chains$K[[1]] <- initialization_values$K
#     chains$G[[1]] <- initialization_values$G
#     chains$z[[1]] <- initialization_values$z 
#     chains$rho[[1]] <- initialization_values$rho 
#     chains$time <- 0
#     chains$sigma <- initialization_values$sigma 
#     chains$theta <- initialization_values$theta
#     chains$weights_a[[1]] <- initialization_values$weights_a
#     chains$weights_d[[1]] <- initialization_values$weights_d
#     # da togliere eventualmente
#     chains$bd_graph_start[[1]] <- initialization_values$graph_start # to check graph G for likelihood of partition
#     chains$total_weights[[1]] <- initialization_values$total_weights # to check graph G for likelihood of partition
#     
#     # initialization of parameters for set_options
#     weights_a <- initialization_values$weights_a
#     weights_d <- initialization_values$weights_d
#     total_weights <- initialization_values$total_weights       
#     total_K <- initialization_values$total_K
#     total_graphs <- initialization_values$total_graphs
#     graph_start <- initialization_values$graph_start
#     
#     seeds = sample(1:9999999, size = (niter + 10))
#     
#     pb = txtProgressBar(min = 2, max = niter, initial = 2, style = 3)  #initialize progress bar
#     
#     for(s in 2:niter) {
#         
#         fit = UpdateParamsGSL(
#             chains$Beta[[s-1]],
#             chains$mu[[s-1]],
#             chains$tau_eps[[s-1]], 
#             chains$K[[s-1]],
#             set_UpdateParamsGSL_list$tbase_base,
#             set_UpdateParamsGSL_list$tbase_data,
#             set_UpdateParamsGSL_list$Sdata,
#             set_UpdateParamsGSL_list$a_tau_eps, 
#             set_UpdateParamsGSL_list$b_tau_eps, 
#             set_UpdateParamsGSL_list$sigma_mu,
#             set_UpdateParamsGSL_list$r,
#             set_UpdateParamsGSL_list$Update_Beta,
#             set_UpdateParamsGSL_list$Update_Mu,
#             set_UpdateParamsGSL_list$Update_Tau,
#             seeds[s]
#         )
#         
#         # Save Beta
#         chains$Beta[[s]] <- fit$Beta
#         
#         # Save mu
#         chains$mu[[s]] <- fit$mu
#         
#         # Save tau
#         chains$tau_eps[[s]] <- fit$tau_eps
#         
#         # Set options for a single iteration of the Gibbs_sampler
#         options = set_options(
#             sigma_prior_0 = chains$sigma[[s-1]],
#             sigma_prior_parameters = list("a"=initialization_values$a_sigma,"b"=initialization_values$b_sigma,
#                                           "c"=initialization_values$c_sigma,"d"=initialization_values$d_sigma),
#             theta_prior_0 = chains$theta[[s-1]],
#             theta_prior_parameters=list("c"=initialization_values$c_theta,"d"=initialization_values$d_theta),
#             rho0=chains$rho[[s-1]],    
#             weights_a0 = weights_a,
#             weights_d0 = weights_d,
#             total_weights0 = total_weights,
#             total_K0 = total_K,
#             total_graphs0 = total_graphs,
#             graph = graph_start,
#             alpha_target = alpha_target,
#             beta_mu = initialization_values$graph_density, 
#             beta_sig2 = initialization_values$beta_sig2, 
#             d = initialization_values$d, 
#             alpha_add = alpha_add, 
#             adaptation_step = adaptation_step,
#             update_sigma_prior = update_sigma_prior,
#             update_theta_prior = update_theta_prior,
#             update_weights = update_weights,
#             update_partition = update_partition,
#             update_graph = update_graph,
#             perform_shuffle = perform_shuffle
#         )
#         
#         # Run an iteration of the Gibbs Sampler
#         res <- Gibbs_sampler(
#             data = t(fit$Beta - fit$mu),
#             niter = 1, nburn = 0, thin = 1,
#             options = options,
#             seed = seed,
#             print = FALSE,
#             s = s # Luciano: Added this for weight adaptivity
#         )
#         
#         z = do.call(rbind, lapply(res$rho, rho_to_z))
#         
#         # Save rho
#         chains$rho[[s]] <- res$rho[[1]]
#         
#         # Save K
#         chains$K[[s]] <- res$K [[1]]
#         
#         # Save G
#         chains$G[[s]] <- res$G [[1]]    
#         
#         # Save z
#         chains$z[[s]] <- z      
#         
#         # Save times for each K
#         chains$time[[s]] <- res$execution_time
#         
#         # Save sigma and theta
#         chains$sigma[[s]] <- res$sigma[[1]]
#         chains$theta[[s]] <- res$theta[[1]]
#         
#         # Save weights
#         chains$weights_a[[s]] <- res$weights_a[[1]]
#         chains$weights_d[[s]] <- res$weights_d[[1]]
#         
#         # Update quantities for the next iteration
#         weights_a <- res$weights_a[[1]]
#         weights_d <- res$weights_d[[1]]
#         total_weights <- res$total_weights
#         total_K <- res$total_K[[1]]
#         total_graphs <- res$total_graphs[[1]]
#         graph_start = res$bdgraph_start
#         
#         # da togliere eventualmente
#         chains$bd_graph_start[[s]] = res$bdgraph_start # to check graph G for likelihood of partition
#         chains$total_weights[[s]] = res$total_weights # to check graph G for likelihood of partition
#         
#         setTxtProgressBar(pb, s)
#     }
#     
#     close(pb)
#     return(chains)
# }


Gibbs_sampler_update = function(
    set_UpdateParamsGSL_list,
    niter, 
    initialization_values, 
    alpha_target,
    alpha_add,
    adaptation_step,
    seed,
    update_sigma_prior,
    update_theta_prior, 
    update_weights,
    update_partition,
    update_graph,
    perform_shuffle,
    algorithm_graph = 'rjmcmc'){
  set.seed(seed)
  # Create a list for chains to save the values of each iteration
  chains <- list(
    Beta = vector("list", length = niter),
    mu = vector("list", length = niter),
    tau_eps = vector("list", length = niter),
    K = vector("list", length = niter),
    G = vector("list", length = niter),
    z = vector("list", length = niter),
    rho = vector("list", length = niter),
    time = vector("list", length = niter), 
    sigma = vector("list", length = niter),
    theta = vector("list", length = niter),
    weights_a = vector("list", length = niter),
    weights_d = vector("list", length = niter),
    # da togliere eventualmente
    bd_graph_start = vector("list", length = niter), 
    total_weights = vector("list", length = niter) 
  ) 
  
  # Initialization of the chains
  chains$Beta[[1]] <- initialization_values$Beta        
  chains$mu[[1]] <- initialization_values$mu
  chains$tau_eps[[1]] <- initialization_values$tau_eps 
  chains$K[[1]] <- initialization_values$K
  chains$G[[1]] <- initialization_values$G
  chains$z[[1]] <- initialization_values$z 
  chains$rho[[1]] <- initialization_values$rho 
  chains$time <- 0
  chains$sigma <- initialization_values$sigma 
  chains$theta <- initialization_values$theta
  chains$weights_a[[1]] <- initialization_values$weights_a
  chains$weights_d[[1]] <- initialization_values$weights_d
  # da togliere eventualmente
  chains$bd_graph_start[[1]] <- initialization_values$graph_start 
  chains$total_weights[[1]] <- initialization_values$total_weights 
  
  # initialization of parameters for set_options
  weights_a <- initialization_values$weights_a
  weights_d <- initialization_values$weights_d
  total_weights <- initialization_values$total_weights       
  total_K <- initialization_values$total_K
  total_graphs <- initialization_values$total_graphs
  graph_start <- initialization_values$graph_start
  
  seeds = sample(1:9999999, size = (niter + 10))
  
  pb = txtProgressBar(min = 2, max = niter, initial = 2, style = 3)  #initialize progress bar
  
  for(s in 2:niter) {
    
    fit = UpdateParamsGSL(
      chains$Beta[[s-1]],
      chains$mu[[s-1]],
      chains$tau_eps[[s-1]], 
      chains$K[[s-1]],
      set_UpdateParamsGSL_list$tbase_base,
      set_UpdateParamsGSL_list$tbase_data,
      set_UpdateParamsGSL_list$Sdata,
      set_UpdateParamsGSL_list$a_tau_eps, 
      set_UpdateParamsGSL_list$b_tau_eps, 
      set_UpdateParamsGSL_list$sigma_mu,
      set_UpdateParamsGSL_list$r,
      set_UpdateParamsGSL_list$Update_Beta,
      set_UpdateParamsGSL_list$Update_Mu,
      set_UpdateParamsGSL_list$Update_Tau,
      seeds[s]
    )
    
    # Save Beta
    chains$Beta[[s]] <- fit$Beta
    
    # Save mu
    chains$mu[[s]] <- fit$mu
    
    # Save tau
    chains$tau_eps[[s]] <- fit$tau_eps
    
    # Set options for a single iteration of the Gibbs_sampler
    options = set_options(
      sigma_prior_0 = chains$sigma[[s-1]],
      sigma_prior_parameters = list("a"=initialization_values$a_sigma,"b"=initialization_values$b_sigma,
                                    "c"=initialization_values$c_sigma,"d"=initialization_values$d_sigma),
      theta_prior_0 = chains$theta[[s-1]],
      theta_prior_parameters=list("c"=initialization_values$c_theta,"d"=initialization_values$d_theta),
      rho0=chains$rho[[s-1]],    
      weights_a0 = weights_a,
      weights_d0 = weights_d,
      total_weights0 = total_weights,
      total_K0 = total_K,
      total_graphs0 = total_graphs,
      graph = graph_start,
      alpha_target = alpha_target,
      beta_mu = initialization_values$graph_density, 
      beta_sig2 = initialization_values$beta_sig2, 
      d = initialization_values$d, 
      alpha_add = alpha_add, 
      adaptation_step = adaptation_step,
      update_sigma_prior = update_sigma_prior,
      update_theta_prior = update_theta_prior,
      update_weights = update_weights,
      update_partition = update_partition,
      update_graph = update_graph,
      perform_shuffle = perform_shuffle
    )
    
    # Run an iteration of the Gibbs Sampler
    res <- Gibbs_sampler(
      data = t(fit$Beta - fit$mu),
      niter = 1, nburn = 0, thin = 1,
      options = options,
      seed = seed,
      print = FALSE,
      s = s,                       # Luciano: Added this for weight adaptivity
      algorithm = algorithm_graph  # Luciano: Added this to choose rjmcmc, rjmcmc.mpl, bdmcmc
    )
    
    z = do.call(rbind, lapply(res$rho, rho_to_z))
    
    # Save rho
    chains$rho[[s]] <- res$rho[[1]]
    
    # Save K
    chains$K[[s]] <- res$K [[1]]
    
    # Save G
    chains$G[[s]] <- res$G [[1]]    
    
    # Save z
    chains$z[[s]] <- z      
    
    # Save times for each K
    chains$time[[s]] <- res$execution_time
    
    # Save sigma and theta
    chains$sigma[[s]] <- res$sigma[[1]]
    chains$theta[[s]] <- res$theta[[1]]
    
    # Save weights
    chains$weights_a[[s]] <- res$weights_a[[1]]
    chains$weights_d[[s]] <- res$weights_d[[1]]
    
    # Update quantities for the next iteration
    weights_a <- res$weights_a[[1]]
    weights_d <- res$weights_d[[1]]
    total_weights <- res$total_weights
    total_K <- res$total_K[[1]]
    total_graphs <- res$total_graphs[[1]]
    graph_start = res$bdgraph_start
    
    # save graphs and respective weights
    chains$bd_graph_start[[s]] = res$bdgraph_start
    chains$total_weights[[s]] = res$total_weights 
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  return(chains)
}


# Added functions for Informed Ordered Partition ---------------------------

# GAMMA THINGS ------------------------------------------------------------

log_A_PY <- function(M_h, p_h, rho_h, theta, sigma) {
  
  # First part: log(n!) - log(K!)
  out <- lfactorial(p_h) - lfactorial(M_h)
  
  # Product over i=1 .. K_h-1
  if (M_h > 1) {
    out <- out + sum(log(theta + sigma * (1:(M_h - 1))))
  }
  
  # # Subtract (K_h - 1) * log((theta+1)_(p_h - 1))
  # if (M_h > 1) {
  #   out <- out - (M_h - 1) * lpochhammer(theta + 1, p_h - 1)
  # }
  
  # Subtract log((theta+1)_(p_h - 1))
  if (M_h > 1) {
    out <- out - lpochhammer(theta + 1, p_h - 1)
  }
  
  # Product over clusters
  for (n_h in rho_h) {
    out <- out + lpochhammer(1 - sigma, n_h - 1) - lfactorial(n_h)
  }
  
  return(out)
}


update_gamma_h_optimized <- function(
    rvec, rvec_current, endpoints,
    modify_endpoint = FALSE,
    gamma,
    eta,
    theta = 2,
    sigma = 0.0001
){
  
  p <- length(rvec)
  gamma_updated <- gamma
  
  # --- update endpoints if needed ---
  if (modify_endpoint) {
    rho_tmp <- r_to_rho(rvec)
    endpoints <- cumsum(rho_tmp)[-length(rho_tmp)]
  }
  
  # --- immediately force incompatible endpoints to gamma=0 ---
  incompatible <- endpoints[ rvec_current[endpoints] == 0 ]
  gamma_updated[incompatible] <- 0
  
  # --- precompute prefix sums once ---
  cs_r <- cumsum(rvec_current)
  
  # helper: fast segment sum
  seg_sum <- function(a, b) {
    if (a > b) return(0L)
    if (a == 1L) return(cs_r[b])
    cs_r[b] - cs_r[a - 1L]
  }
  
  # --- helper: fast rho computation from changepoints ---
  fast_rho <- function(a, b) {
    # a..b is a small contiguous region
    cp <- which(rvec_current[a:b] == 1)
    if (length(cp) == 0L) return( integer(0) )
    diff(c(0L, cp))
  }
  
  # diagnostics
  diag <- vector("list", length(endpoints))
  names(diag) <- endpoints
  
  # MAIN LOOP
  for (node in endpoints) {
    
    # skip nodes forced to 0
    if (gamma_updated[node] == 0 && node %in% incompatible)
      next
    
    # Case γ_j = 1  (numerator)
    gamma_mod1 <- gamma
    gamma_mod1[node] <- 1L
    gpos1 <- which(gamma_mod1 == 1L)
    
    # previous/next gamma positions
    g_before <- gpos1[gpos1 < node]
    g_after  <- gpos1[gpos1 > node]
    
    # h_minus range
    if (length(g_before) == 0L) {
      a_minus <- 1L
    } else {
      a_minus <- max(g_before) + 1L
    }
    b_minus <- node
    
    # h_plus range
    if (length(g_after) == 0L) stop("No gamma after node in numerator")
    a_plus <- node + 1L
    b_plus <- g_after[1L]
    
    # statistics
    M_minus <- seg_sum(a_minus, b_minus)
    M_plus  <- seg_sum(a_plus, b_plus)
    
    p_minus <- b_minus - a_minus + 1L
    p_plus  <- b_plus  - a_plus  + 1L
    
    rho_minus <- fast_rho(a_minus, b_minus)
    rho_plus  <- fast_rho(a_plus, b_plus)
    
    # Case γ_j = 0 (denominator)
    gamma_mod0 <- gamma
    gamma_mod0[node] <- 0L
    gpos0 <- which(gamma_mod0 == 1L)
    
    g_before0 <- gpos0[gpos0 < node]
    g_after0  <- gpos0[gpos0 > node]
    
    if (length(g_after0) == 0L) stop("No gamma after node in denominator")
    
    # low part
    if (length(g_before0) == 0L) {
      a_low <- 1L
    } else {
      a_low <- max(g_before0) + 1L
    }
    b_low <- node
    
    # high part
    a_high <- node + 1L
    b_high <- g_after0[1L]
    
    M_u <- seg_sum(a_low, b_high)
    p_u <- b_high - a_low + 1L
    rho_u <- fast_rho(a_low, b_high)
    
    # Compute log terms
    log_A_minus <- log_A_PY(M_minus, p_minus, rho_minus, theta, sigma)
    log_A_plus  <- log_A_PY(M_plus,  p_plus,  rho_plus,  theta, sigma)
    log_B_u     <- log_A_PY(M_u,     p_u,     rho_u,     theta, sigma)
    
    lnN <- log(eta)     + log_A_minus + log_A_plus
    lnD <- log(1 - eta) + log_B_u
    
    P1 <- plogis(lnN - lnD)
    
    # Gibbs update
    gamma_updated[node] <- rbinom(1L, 1L, P1)
    
    diag[[as.character(node)]] <- data.frame(
      node = node,
      lnN = lnN,
      lnD = lnD,
      num = exp(lnN),
      den = exp(lnN) + exp(lnD),
      P1 = P1,
      log_eta = log(eta),
      log_A_minus = log_A_minus,
      log_A_plus  = log_A_plus,
      log_1_min_eta = log(1 - eta),
      log_B_u     = log_B_u,
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    gamma_updated = gamma_updated,
    diagnostics = diag
  ))
}


# PARTITION THINGS --------------------------------------------------------
# split into subgroups according to gamma

split_by_gamma <- function(rho, gamma) {
  cum_sizes <- cumsum(rho)
  ends_gamma <- which(gamma == 1)
  ends_blocks <- findInterval(ends_gamma, cum_sizes)
  starts_blocks <- c(1, ends_blocks[-length(ends_blocks)] + 1)
  mapply(function(a, b) rho[a:b], starts_blocks, ends_blocks, SIMPLIFY = FALSE)
}


split_weights_by_gamma <- function(w, gamma) {
  if(length(w) != length(gamma)) stop("wa and gamma must have the same length")
  # Find the positions of gamma == 1
  cp <- which(gamma == 1)
  # Start positions: first element and one after each previous changepoint
  starts <- c(1, cp[-length(cp)] + 1)
  ends <- cp
  # Split wa into sublists
  w_list <- mapply(function(a, b) w[a:b], starts, ends, SIMPLIFY = FALSE)
  return(w_list)
}



proposal_ratio_h = function(rho_h,  
                            alpha_add,
                            wa_h,
                            wd_h,
                            choose_add,
                            h, H, p) {
  # number of innergroups in subgroup h and nodes in subgroup h
  M_h = length(rho_h)
  p_h = sum(rho_h)
  
  # indexes of the changepoints in sub-group h
  cp_indexes_h <- get_group_indexes(rho_h)
  
  # if((h==H) & (length(cp_indexes_h) == 1) & (cp_indexes_h[length(cp_indexes_h)]==p_h)){ 
  #   # if we are in the last subgroup, H, and only one node we do not update
  #   # cp_indexes_h <- cp_indexes_h[-length(cp_indexes_h)]
  #   # candidate = rho_h
  #   return(list("ratio" = NULL, "candidate" = 'last_node'))
  #   }
  
  # if(h==H && length(cp_indexes_h) == 0){ # if only one node in last sub-group
  #   # then let it a changepoint
  #   # candidate = rho_h
  #   return(list("ratio" = NULL, "candidate" = 'last_node'))
  # }
  
  # not all points can be selected for add move --> assign probability zero to those who cannot be
  wa_h_available = wa_h
  wa_h_available[cp_indexes_h] = 0
  wa_h_available_sum = sum(wa_h_available)
  
  # not all points can be selected for delete move --> assign probability zero to those who cannot be
  wd_h_available = wd_h
  wd_h_available[-cp_indexes_h] = 0
  wd_h_available_sum = sum(wd_h_available)
  
  if (choose_add) {
    draw_weights = wa_h_available
  } else {
    draw_weights = wd_h_available
  }
  
  if(
    ((cp_indexes_h[length(cp_indexes_h)]==p) & (length(draw_weights) != p)) | #rho=p and draw_weights=39
    ((h==H) & (length(draw_weights) != p_h))  # more sub-groups and at the last we choosen to delete a cpn
  ){
    # draw the candidate among the first 1:(p-1)
    candidate = sample(1:(p_h-1), 1, prob = draw_weights)
  }else{
    candidate = sample(1:p_h, 1, prob = draw_weights) 
  }
  
  
  if (choose_add && M_h == 1) {
    # case in which you choose to propose an add move
    # (with just 1 group) that may or may not be accepted
    ratio = ((1-alpha_add)/alpha_add) * (wa_h_available_sum / wa_h[candidate])
    return(list("ratio" = ratio, "candidate" = candidate#, 
                #"wa_h_available" = wa_h_available,
                #"wd_h_available" = wd_h_available
    ))
  }
  
  if (!choose_add && M_h == p_h) {
    # case in which you choose to propose a delete move
    # (with every point being a group) that may or may not be accepted
    ratio = (alpha_add /(1-alpha_add)) * (wd_h_available_sum / wd_h[candidate])
    return(list("ratio" = ratio, "candidate" = candidate#, 
                #"wa_h_available" = wa_h_available,
                #"wd_h_available" = wd_h_available
    ))
  }
  
  # only the general cases remain
  if (choose_add) {
    ratio = (1 - alpha_add) / alpha_add *
      (wa_h_available_sum / wa_h[candidate]) *
      ( wd_h[candidate] / (wd_h[candidate] + wd_h_available_sum) )
  } else {
    ratio = alpha_add / (1 - alpha_add) *
      (wd_h_available_sum / wd_h[candidate]) *
      ( wa_h[candidate] / (wa_h[candidate] + wa_h_available_sum) )
  }
  
  return(list("ratio" = ratio, "candidate" = candidate#, 
              #"wa_h_available" = wa_h_available,
              #"wd_h_available" = wd_h_available
  ))
}


split_partition_h = function(candidate_index, rho_h) {
  
  # number of groups
  M_h = length(rho_h)
  new_rho_h = rep(NA, M_h + 1)
  
  group_indexes = get_group_indexes(rho_h)
  found = FALSE
  
  for (i in 1:M_h) {
    
    # update the partition in the general case
    # (either I have already split the group or not, just the index changes)
    if (!found) {
      new_rho_h[i] = rho_h[i]
    } else {
      new_rho_h[i + 1] = rho_h[i]
    }
    
    if (!found && group_indexes[i] > candidate_index) {
      # just passed the element index - I am in the group to be split
      
      # index of the element minus the cumulative
      # number of elements in the previous groups only if i!=1
      new_rho_h[i] = candidate_index - (i != 1) * group_indexes[i - 1 * (i != 1)]
      # dimension of the original group minus the elements moved to new_rho[i]
      new_rho_h[i + 1] = rho_h[i] - new_rho_h[i]
      
      # save the index of the group that has changed
      j = i
      
      found = TRUE
    }
  }
  return(list("new_rho_h" = new_rho_h, "changed_group_index" = j))
}

merge_partition_h = function(candidate_index, rho_h) {
  
  # number of groups
  M_h= length(rho_h)
  new_rho_h = rep(NA, M_h - 1)
  
  group_indexes = get_group_indexes(rho_h)
  found = FALSE
  
  # Valid merge points are group_indexes[1:(M_h-1)]
  valid_merge_points <- group_indexes[1:(M_h - 1)]
  if (!(candidate_index %in% valid_merge_points)) {
    return(list("new_rho_h" = rho_h,
                "changed_group_index" = NA,
                invalid_move = TRUE))
  }
  
  for (i in 1:(M_h - 1)) {
    
    # update the partition in the general case
    # either I have already merged the group or not, just the index changes
    if (!found) {
      new_rho_h[i] = rho_h[i]
    } else {
      new_rho_h[i] = rho_h[i + 1]
    }
    
    if (!found && group_indexes[i] == candidate_index) {
      # I am at the changepoint between the two groups to be merged
      
      # index of the element minus the cumulative
      # number of elements in the previous groups
      new_rho_h[i] = rho_h[i] + rho_h[i + 1]
      
      # save the index of the group that has changed
      j = i
      
      found = TRUE
    }
  }
  return(list("new_rho_h" = new_rho_h,
              "changed_group_index" = j,
              invalid_move = FALSE))
}


# rho_h_current = rho_h # needed to see how function works
log_prior_ratio_h = function(theta_prior,
                             sigma_prior,
                             rho_h_current,
                             rho_h_proposed,
                             choose_add,
                             changed_group_index)
{ 
  # differentiate delete/merge case
  if (!choose_add) {
    # swap rhos 'cause we're lazy
    temp = rho_h_current
    rho_h_current = rho_h_proposed
    rho_h_proposed = temp
  }
  
  # number of groups in the current partition
  M_h = length(rho_h_current)
  
  s_h = changed_group_index
  #s = get_index_changed_group(rho_h_current,rho_h_proposed)
  
  # compute the prior ratio
  log_ratio = - log(M_h+1) + log(theta_prior + M_h * sigma_prior)
  + lpochhammer(1 - sigma_prior, rho_h_proposed[s_h] - 1)
  + lpochhammer(1 - sigma_prior, rho_h_proposed[s_h + 1] - 1)
  - lpochhammer(1 - sigma_prior, rho_h_proposed[s_h] + rho_h_proposed[s_h + 1] - 1)
  + lfactorial(rho_h_proposed[s_h] + rho_h_proposed[s_h + 1])
  - lfactorial(rho_h_proposed[s_h])
  - lfactorial(rho_h_proposed[s_h + 1])
  
  
  # in the delete/merge case we have to invert everything
  if (!choose_add) {
    log_ratio = -log_ratio
  }
  
  return(log_ratio)
}



# # log lik ratio
# a1=estimate_Beta_params(mu=0.8,var=0.01)[[1]]
# a2=estimate_Beta_params(mu=0.8,var=0.01)[[2]]
# # G_last = tail(chains$G, n=1)[[1]]
# # G_list = BFDR_selection(G_last, tol = seq(0.1, 1, by = 0.001))
# # G = G_list$best_truncated_graph
# G = G_est
# rho_h_current = rho_h
# rho_current = rho_updated # # this is rho at iter t-1, so it is rho_current at t
# rho_proposed <- rho_current_list
# rho_proposed[[h]] <- rho_h_proposed
# rho_proposed = unlist(rho_proposed)

log_likelihood_ratio_h <- function(alpha_add,
                                   wa_h,
                                   wd_h,
                                   G,
                                   #rho_h_current,## to keep??
                                   rho_current,
                                   rho_current_list,
                                   #rho_h_proposed, ## to keep??
                                   rho_proposed,
                                   choose_add,
                                   a1,
                                   a2,
                                   changed_group_index,
                                   h, H) {
  # differentiate delete/merge case
  if (!choose_add) {
    # swap rhos 'cause we're lazy
    temp = rho_current
    rho_current = rho_proposed
    rho_proposed = temp
  }
  
  # number of groups
  M = length(rho_current)
  
  # wrap the general fB into this version so prior beta params are specified once 
  fB = function(group1, group2, S, S_star) {
    return(fB_general(
      group1,
      group2,
      S,
      S_star,
      alpha = a1,
      beta = a2,
      log = TRUE
    ))
  }
  
  S_current = get_S_from_G_rho(G, rho_current)
  S_star_current = get_S_star_from_S_and_rho(S_current, rho_current)
  
  S_proposed = get_S_from_G_rho(G, rho_proposed)
  S_star_proposed = get_S_star_from_S_and_rho(S_proposed, rho_proposed)
  
  # index of group that changed in current rho within sub-group h
  s_h = changed_group_index 
  # index of group that changed in whole current partition 
  s_overall <- sum(lengths(rho_current_list)[seq_len(h - 1)]) + s_h
  
  #log of ( 1/B(a1,a2) )^(K+1) in the paper
  log_ratio = -(M + 1) * fB_zero(alpha = a1, beta = a2)
  
  # for rows < s_overall, col = s_overall and s_overall + 1 
  for (l in 1:(s_overall - 1)) {
    if (l > (s_overall - 1)) break; # avoiding strange things when s_overall and M_h are too small
    # numerator term
    log_ratio = log_ratio + fB(l, s_overall, S_proposed, S_star_proposed)
    log_ratio = log_ratio + fB(l, s_overall + 1, S_proposed, S_star_proposed)
    # denominator term
    log_ratio = log_ratio - fB(l, s_overall, S_current, S_star_current)
  }
  
  # numerator, for rows = s_overall and s_overall + 1, col >= s_overall + 2 
  for (m in (s_overall + 2):(M + 1)) {
    if (m > (M + 1)) break; # avoiding strange things when s_overall and M are too small
    log_ratio = log_ratio + fB(s_overall, m, S_proposed, S_star_proposed) +
      fB(s_overall + 1, m, S_proposed, S_star_proposed)
  }
  
  # denominator, for rows = s_overall, col >= s_overall + 1 (or col - 1 >= s_overall + 2) 
  for (m in (s_overall + 1):M) {
    if (m > M) break; # avoiding strange things when s_overall and M_h are too small
    log_ratio = log_ratio - fB(s_overall, m, S_current, S_star_current)
  }
  
  # numerator, for rows = s_overall and s_overall + 1, col = s_overall and s_overall + 1
  log_ratio = log_ratio + fB(s_overall, s_overall + 1, S_proposed, S_star_proposed) +
    fB(s_overall, s_overall, S_proposed, S_star_proposed) +
    fB(s_overall + 1, s_overall + 1, S_proposed, S_star_proposed)
  
  # denominator, for rows = s_overall, col = s_overall 
  log_ratio = log_ratio - fB(s_overall, s_overall, S_current, S_star_current)
  
  # in the delete/merge case we have to invert everything
  if (!choose_add) {
    log_ratio = -log_ratio
  }
  
  return(log_ratio)
  
}


update_partition_h = function(rho_current_list, # defined outside for loop on h=1...H in Gibbs_sampler_h 
                              rho_h_current,
                              alpha_add,
                              wa_h,
                              wd_h,
                              theta_prior,
                              sigma_prior,
                              G,
                              beta_params, 
                              h, H, p){
  
  rho_current = unlist(rho_current_list) # this remains as it is in the for loop h=1...H
  
  unifsample = runif(n = 1)
  choose_add = unifsample < alpha_add
  M_h = length(rho_h_current)  # number of groups, in paper K_h
  p_h = sum(rho_h_current)     # number of nodes into subgroup, in paper tilde_n_h
  # stop and return current if only one node
  if((M_h==p_h)&(p_h==1)){
    return(list(
      "rho_updated" = rho_h_current,
      "accepted" = as.numeric(FALSE),
      "choose_add" = choose_add,
      "candidate" = M_h,
      "alpha_accept" = 0
    )) 
  }
  # change if not feasibility of choose_add
  if ((!choose_add && M_h == 1) || (choose_add && M_h == p_h)) { 
    choose_add = !choose_add
  }
  choose_add
  proposal_list_h = proposal_ratio_h(rho_h = rho_h_current, 
                                     alpha_add = alpha_add, 
                                     wa_h = wa_h, wd_h = wd_h, 
                                     choose_add = choose_add,
                                     h = h, H = H, p=p)
  proposal_list_h
  
  # if it is the last node than we do not modify anything 
  if(proposal_list_h$candidate == "last_node"){
    return(list(
      "rho_updated" = rho_h_current,
      "accepted" = as.numeric(FALSE),
      "choose_add" = choose_add,
      "candidate" = rho_h_current,
      "alpha_accept" = 0
    ))
  }
  log_proposal_ratio_h = log(proposal_list_h$ratio)
  candidate_h = proposal_list_h$candidate # this is the node
  # wa_h_available = proposal_list_h$wa_h_available
  # wd_h_available = proposal_list_h$wd_h_available
  if (choose_add) {
    list_output_modify_partition_h = split_partition_h(candidate_index = candidate_h, 
                                                       rho_h = rho_h_current#,
                                                       #wa_h_available = wa_h_available
    )
  } else {
    list_output_modify_partition_h = merge_partition_h(candidate_index = candidate_h,
                                                       rho_h = rho_h_current#,
                                                       #wd_h_available = wd_h_available
    )
  }
  # if we have an invalid proposal wrt the subgroup h (i.e. if propose for 
  # merging the last node in h) then we reject and keep old partition
  if (!is.null(list_output_modify_partition_h$invalid_move) &&
      list_output_modify_partition_h$invalid_move) {
    # Immediately reject the proposal
    accepted = FALSE
    return(list(
      "rho_updated" = rho_h_current,
      "accepted" = as.numeric(accepted),
      "choose_add" = choose_add,
      "candidate" = candidate_h,
      "alpha_accept" = 0
    ))
  }
  
  # Otherwise we keep to see the proposal
  rho_h_proposed  = list_output_modify_partition_h$new_rho_h
  
  # changing the sub-group h from the current partition to the proposed partition
  rho_proposed <- rho_current_list
  rho_proposed[[h]] <- rho_h_proposed
  rho_proposed = unlist(rho_proposed)
  
  changed_group_index = list_output_modify_partition_h$changed_group_index # this is s_h in following functions
  
  l_prior_ratio_h = log_prior_ratio_h(
    theta_prior,
    sigma_prior,
    rho_h_current = rho_h_current,
    rho_h_proposed = rho_h_proposed,
    choose_add,
    changed_group_index
  )
  
  l_likelihood_ratio_h = log_likelihood_ratio_h(
    alpha_add,
    wa_h = wa_h,
    wd_h = wd_h,
    G = G, # this is the graph in the update_partition_h parameter
    #rho_h_current = rho_h_current,
    rho_current,
    rho_current_list,
    #rho_h_proposed = rho_h_proposed,
    rho_proposed,
    choose_add,
    a1 = beta_params$alpha,
    a2 = beta_params$beta,
    changed_group_index,
    h = h, H = H
  )
  
  alpha_accept <- min(1, exp(l_likelihood_ratio_h +
                               l_prior_ratio_h +
                               log_proposal_ratio_h))
  if (runif(n = 1) < alpha_accept) {
    accepted = TRUE
    rho_h_updated = rho_h_proposed
  } else {
    accepted = FALSE
    rho_h_updated = rho_h_current
  }
  # accepted
  # rho_h_updated
  return(
    list(
      "rho_updated" = rho_h_updated,
      "accepted" = as.numeric(accepted),
      "choose_add" = choose_add,
      "candidate" = candidate_h,
      "alpha_accept" = alpha_accept
    )
  )
}



set_initialization_h <- function(
    Beta,
    mu,
    tau_eps,
    K,
    G,
    z,
    rho,
    a_sigma,
    b_sigma,
    c_sigma,
    d_sigma,
    c_theta,
    d_theta,
    sigma,
    theta,
    weights_a0,
    weights_d0,
    total_weights,
    total_K,
    total_graphs,
    graph_start,
    graph_density,
    beta_sig2,
    d,
    gamma
){
  
  initialization_values = list(
    'Beta'            = Beta,
    'mu'              = mu,
    'tau_eps'         = tau_eps,
    'K'               = K,
    'G'               = G,
    'z'               = z,
    'rho'             = rho,
    'a_sigma'         = a_sigma,
    'b_sigma'         = b_sigma,
    'c_sigma'         = c_sigma,
    'd_sigma'         = d_sigma,
    'c_theta'         = c_theta,
    'd_theta'         = d_theta,
    'sigma'           = sigma,
    'theta'           = theta,
    'weights_a'       = weights_a0,
    'weights_d'       = weights_d0,
    'total_weights'   = total_weights,
    'total_K'         = total_K,
    'total_graphs'    = total_graphs,
    'graph_start'     = graph_start,
    'graph_density'   = graph_density,
    'beta_sig2'       = beta_sig2,
    'd'               = d,
    'gamma'           = gamma
  )
  
  return(initialization_values)
}



set_options_h <- function(sigma_prior_0,
                          sigma_prior_parameters,
                          theta_prior_0,
                          theta_prior_parameters,
                          rho,
                          weights_a0,
                          weights_d0,
                          total_weights0, 
                          total_K0, 
                          total_graphs0, 
                          graph,   
                          alpha_target,
                          beta_mu,
                          beta_sig2,
                          d=3,
                          alpha_add=0.5,
                          adaptation_step,
                          update_sigma_prior=TRUE,
                          update_theta_prior=TRUE,
                          update_weights=TRUE,
                          update_partition=TRUE,
                          update_graph=TRUE,
                          perform_shuffle=TRUE,
                          # informed partition
                          gamma,
                          update_gamma = T,
                          compute_partition_update_info = FALSE,
                          eta,
                          sample_eta) {
  
  options = list(
    "sigma_prior_0"          = sigma_prior_0,
    "sigma_prior_parameters" = sigma_prior_parameters,
    "theta_prior_0"          = theta_prior_0,
    "theta_prior_parameters" = theta_prior_parameters,
    "rho"                    = rho,
    "weights_a0"             = weights_a0,
    "weights_d0"             = weights_d0,
    "total_weights0"         = total_weights0,
    "total_K0"               = total_K0,
    "total_graphs0"          = total_graphs0,
    "graph"                  = graph,
    "alpha_target"           = alpha_target,
    "beta_mu"                = beta_mu,
    "beta_sig2"              = beta_sig2,
    "d"                      = d,
    "alpha_add"              = alpha_add,
    "adaptation_step"        = adaptation_step,
    "update_sigma_prior"     = update_sigma_prior,
    "update_theta_prior"     = update_theta_prior,
    "update_weights"         = update_weights,
    "update_partition"       = update_partition,
    "update_graph"           = update_graph,
    "perform_shuffle"        = perform_shuffle,
    
    'gamma'                  = gamma,
    'update_gamma'           = update_gamma,
    'compute_partition_update_info'  = compute_partition_update_info,
    'eta'                    = eta,
    'sample_eta'             = sample_eta
  )
  return(options)
}

# NEW GIBBS INNER ---------------------------------------------------------

Gibbs_sampler_h = function(data,
                           niter,
                           nburn,
                           thin,
                           options,
                           seed=1234,
                           print=TRUE,
                           s,
                           rho_0,
                           algorithm)
{
  n = nrow(data) # number of observations
  p = ncol(data) # number of nodes
  n_total_iter = nburn + niter * thin # total iterations to be made
  
  # dynamic parameters
  sigma_prior            = options$sigma_prior_0
  theta_prior            = options$theta_prior_0
  rho                    = options$rho
  weights_a              = options$weights_a0
  weights_d              = options$weights_d0
  total_weights          = options$total_weights0     # added
  total_K                = options$total_K0           # added
  total_graphs           = options$total_graphs0      # added
  graph                  = options$graph              # added
  last_S                 = NULL
  adaptation_step        = options$adaptation_step
  sigma_prior_parameters = options$sigma_prior_parameters
  theta_prior_parameters = options$theta_prior_parameters
  
  # informed partition parameters
  gamma = options$gamma
  rvec = c(rho_to_r(rho_0),1) # initial expert partition in cpn format
  internal_endpoints <- cumsum(rho_0)[-length(rho_0)] # expert cpns
  compute_partition_update_info <- options$compute_partition_update_info # for debug df 
  eta = options$eta
  # constant parameters
  alpha_add    = options$alpha_add
  alpha_target = options$alpha_target
  
  # parameter for the Wishart
  d = options$d
  
  # parameters for the Beta
  beta_mu   = options$beta_mu
  beta_sig2 = options$beta_sig2
  
  
  # checks
  if(sum(rho) != p)
    stop("The partition rho must sum to the number of variables p")
  if(d < 3)
    stop("The Wishart's d parameter must be greater or equal than 3")
  if(!(beta_mu > 0 && beta_mu < 1))
    stop("The mean of the Beta must be between 0 and 1")
  if(!(beta_sig2 > 0 && beta_sig2 < beta_mu*(1-beta_mu)))
    stop("The variance of the Beta must be between 0 and beta_mu*(1-beta_mu)")
  if(length(weights_a) != (p-1) || length(weights_d) != (p-1))
    stop("The number of elements in the weights vectors must be equal to p-1")
  if(!(adaptation_step > 0))
    stop("The adapation step h must be positive")
  if(!(alpha_add > 0 && alpha_add < 1))
    stop("The probability of choosing an add move alpha_add must be between 0 and 1")
  if(!(alpha_target > 0 && alpha_target < 1))
    stop("The target acceptance rate of the Metropolis-Hastings alpha_target must be between 0 and 1")
  
  # t_over_p = n_total_iter / p
  t_over_p = s / p # changed t_over_p from upper line to this line for adaptivity
  
  beta_params = estimate_Beta_params(beta_mu, beta_sig2)
  
  # define structure to save sampled values
  save_res = list(
    G = vector("list", length = niter),
    K = vector("list", length = niter),
    rho = vector("list", length = niter),
    accepted = vector("numeric", length = niter),
    S = vector("list", length = niter),
    sigma = vector("numeric", length = niter),
    theta = vector("numeric", length = niter),
    weights_a = vector("list", length = niter),
    weights_d = vector("list", length = niter),
    bdgraph_start = NULL,
    total_graphs = vector("list", length = niter),
    total_K = vector("list", length = niter),
    total_weights = 0,
    
    gamma = vector("list", length = niter),
    partition_update_info = vector("list", niter),
    eta = vector("numeric", length = niter)
  )
  
  # initialize iteration counter
  it_saved = 0
  
  # initialize progress bar
  if(print){
    pb = txtProgressBar(min=1, max=n_total_iter, initial=1, style=3)
  }
  
  # save start time for measuring execution time
  start_time = Sys.time()
  
  partition_update_info_df <- NULL
  if (compute_partition_update_info) {
    partition_update_info_df <- data.frame(
      iter = integer(0),
      h = integer(0),
      accepted = numeric(0),
      choose_add = logical(0),
      candidate = integer(0)
    )
  }
  
  # start the simulation
  for(iter in 1:n_total_iter){
    
    if(is.null(graph)){ 
      # we run a single iteration of BDgraph with iter = 1 and burnin = 0
      # and g.start = "empty"
      output = post_graph_sampling(
        data,
        rho,
        n,
        method = "ggm",
        algorithm = algorithm,
        iter = 1,
        burnin = 0,
        not.cont = NULL,
        g.prior = 0.5,
        df.prior = d,
        CCG_D = NULL,
        g.start = "empty",
        jump = NULL,
        save = TRUE,
        print = 1,
        cores = NULL,
        threshold = 1e-8,
        beta_params = beta_params
      )
    }
    
    else {
      # we run a single iteration of BDgraph with iter = 1 and burnin = 0
      # and g.start = graph
      output = post_graph_sampling(
        data,
        rho,
        n,
        method = "ggm",
        algorithm = algorithm,
        iter = 1,
        burnin = 0,
        not.cont = NULL,
        g.prior = 0.5,
        df.prior = d,
        CCG_D = NULL,
        g.start = graph,
        jump = NULL,
        save = TRUE,
        print = 1,
        cores = NULL,
        threshold = 1e-8,
        beta_params = beta_params
      )
    }
    
    # update graph
    if (options$update_graph){
      
      # extract precision matrix K
      last_K = output$last_K
      
      if(algorithm == 'bdmcmc'){ 
        # update total_weights
        total_weights = total_weights + output$all_weights
        # update total_graphs taking into consideration the weights
        total_graphs = total_graphs + output$last_graph * output$all_weights
        total_K = total_K + output$last_K * output$all_weights
      }else{
        # in this case we are using rjmcmc hence every graph and precision have weight 1 
        # if doing 1 iteration
        total_weights = total_weights + output$all_weights 
        total_graphs = total_graphs + output$last_graph
        total_K = total_K + output$last_K
      }
    }else{
      # extract precision matrix K even if not updating graph
      last_K = output$last_K  
    }
    
    # extract adjacency matrix G and update for the next iteration
    graph = output$last_graph
    
    # NEW options$update_gamma ----------------------------------
    if (!options$update_gamma) {
      gamma_updated <- gamma
    }
    
    if(options$update_gamma){
      if(options$sample_eta){
        eta_sampled = runif(1,min = 0, max = 0.99)
        #eta_sampled = rbeta(1,0.01,1)
        #eta_sampled = rbeta(1,0.1,1)
      } else{
        eta_sampled = eta
      }
      rvec_current = c(rho_to_r(rho),1)
      gamma_updated_list <- update_gamma_h_optimized(rvec = rvec, 
                                                     rvec_current = rvec_current, 
                                                     endpoints = internal_endpoints, 
                                                     modify_endpoint = F, 
                                                     gamma = gamma, 
                                                     eta = eta_sampled, 
                                                     theta = theta_prior, 
                                                     sigma = sigma_prior)
      gamma_updated = gamma_updated_list$gamma_updated
    }
    
    
    # NEW options$update_partition -------------------------------
    if(options$update_partition){
      # adding 1 because I coded everything of partition update with last node included
      weights_a = c(weights_a,1) 
      weights_d = c(weights_d,1) 
      
      rho_current_list = split_by_gamma(rho, gamma_updated)
      #print(paste0('rho current list: ', rho_current_list))
      wa_list  = split_weights_by_gamma(weights_a, gamma_updated)
      wd_list  = split_weights_by_gamma(weights_d, gamma_updated)
      
      H = sum(gamma_updated) # number of subgroups
      if(!H == length(rho_current_list)){
        stop('H is not of the same length of the current rho list')
      }
      
      # Drop final position in wa_H, wd_H 
      wa_list[[H]] = wa_list[[H]][-length(wa_list[[H]])]
      wd_list[[H]] = wd_list[[H]][-length(wd_list[[H]])]
      
      # Iteration over sub-groups h=1,..,H given by gamma
      rho_list_updated <- vector("list", H)
      wa_h_updated <- vector("list", H)
      wd_h_updated <- vector("list", H)
      
      iter_partition_info <- vector("list", H)
      
      for(h in 1:H){
        
        rho_h = rho_current_list[[h]]
        wa_h  = wa_list[[h]]
        wd_h  = wd_list[[h]]
        
        uph = update_partition_h(
          rho_current_list = rho_current_list,
          rho_h, alpha_add, wa_h, wd_h,
          theta_prior, sigma_prior,
          graph, 
          beta_params, h = h, H = H, p=p)
        
        alpha_accept = uph$alpha_accept
        
        # update of weights
        if(options$update_weights){
          #if(T)
          
          if(uph$accepted){ 
            if(uph$choose_add){
              wa_h = update_weight(
                weights = wa_h,
                uph$candidate,
                adaptation_step,
                t_over_p,
                alpha_accept,
                alpha_target
              )
            } else {
              wd_h = update_weight(
                weights = wd_h,
                uph$candidate,
                adaptation_step,
                t_over_p,
                alpha_accept,
                alpha_target
              )
            }
          }  
        }# end weights update
        
        rho_list_updated[[h]] = uph$rho_updated
        wa_h_updated[[h]] = wa_h 
        wd_h_updated[[h]] = wd_h
        
        iter_partition_info[[h]] <- list(
          accepted = uph$accepted,
          choose_add = uph$choose_add,
          candidate = uph$candidate
        )
        
      } # end iteration over sub-groups h=1...H
      
      weights_a <- unlist(wa_h_updated)
      weights_d <- unlist(wd_h_updated)
      # weights_a = c(weights_a,1) # don't do this
      # weights_d = c(weights_d,1) # don't do this
      # it is here because of the way I coded the update partition, separately 
      # from this Gibbs, but here not necessary to add the last 1, just add it 
      # at the start of the update partition
      rho <- unlist(rho_list_updated)
    } # end update partition

  
    if(options$perform_shuffle){
      rho = shuffle_partition(rho, graph, sigma_prior, beta_params$alpha, beta_params$beta)
    }

    
    if(options$update_sigma_prior){
      candidate <- runif(1,max(0,-theta_prior),1)
      alpha_MH <- full_conditional_sigma(candidate,
                                         theta_prior,
                                         rho,
                                         sigma_prior_parameters$a,
                                         sigma_prior_parameters$b,
                                         sigma_prior_parameters$c,
                                         sigma_prior_parameters$d) -
        full_conditional_sigma(sigma_prior,
                               theta_prior,
                               rho,
                               sigma_prior_parameters$a,
                               sigma_prior_parameters$b,
                               sigma_prior_parameters$c,
                               sigma_prior_parameters$d)
      
      if(log(runif(1)) <= min(alpha_MH,log(1))){
        sigma_prior = candidate
      }else{
        sigma_prior = sigma_prior
      }
    }
    
    if(options$update_theta_prior) {
      theta_prior = full_conditional_theta(
        theta_prior_parameters$c,
        theta_prior_parameters$d,
        theta_prior,
        length(rho),
        p,
        sigma_prior
      )
    }
    
    if(options$update_graph){
      last_S = get_S_from_G_rho(graph,rho)
    }
    
    # save results only on thin iterations
    # (i.e. only save multiples of thin)
    if(iter > nburn && (iter - nburn) %% thin == 0) {
      it_saved = it_saved + 1
      if(options$update_graph){
        # cumulative precision matrix K and probability of inclusion links
        save_res$K[[it_saved]] = total_K / total_weights
        save_res$G[[it_saved]] = total_graphs / total_weights
      }
      else{
        save_res$K[[it_saved]] = total_K
        save_res$G[[it_saved]] = total_graphs
      }
      
      save_res$total_weights <- total_weights
      save_res$total_K[[it_saved]] <- total_K
      save_res$total_graphs[[it_saved]] <- total_graphs
      
      # save bdgraph object
      save_res$bdgraph_start = output$last_graph
      
      save_res$rho[[it_saved]] = rho
      save_res$weights_d[[it_saved]] <- weights_d
      save_res$weights_a[[it_saved]] <- weights_a
      
      save_res$gamma[[it_saved]] <- gamma_updated # informed partition
      
      if(options$sample_eta){
        save_res$eta[[it_saved]] = eta_sampled
      }else{
        save_res$eta[[it_saved]] = eta_sampled
      }
      
      # partition_update_info[[iter]] <- iter_partition_info # informed partition
      # 
      # partition_update_info_df <- do.call(
      #   rbind,
      #   lapply(seq_along(partition_update_info), function(iter) {
      #     iter_info <- partition_update_info[[iter]]
      #     do.call(
      #       rbind,
      #       lapply(seq_along(iter_info), function(h) {
      #         data.frame(
      #           iter = iter,
      #           h = h,
      #           accepted = iter_info[[h]]$accepted,
      #           choose_add = iter_info[[h]]$choose_add,
      #           candidate = iter_info[[h]]$candidate
      #         )
      #       })
      #     )
      #   })
      # )
      
      if (compute_partition_update_info) {
        iter_info_df <- data.frame(
          iter = s,
          h = seq_along(iter_partition_info),
          accepted = vapply(iter_partition_info, `[[`, 0, "accepted"),
          choose_add = vapply(iter_partition_info, `[[`, TRUE, "choose_add"),
          candidate = vapply(iter_partition_info, `[[`, NA_integer_, "candidate")
        )
        partition_update_info_df <- rbind(partition_update_info_df, iter_info_df)
        save_res$partition_update_info <- partition_update_info_df
        save_res$accepted <- sum(partition_update_info_df$accepted)
      } else {
        save_res$partition_update_info <- NULL
        save_res$accepted = sum(vapply(iter_partition_info, `[[`, 0, "accepted"))
      }
      
      save_res$sigma[[it_saved]] = sigma_prior
      save_res$theta[[it_saved]] = theta_prior
      save_res$S[[it_saved]] = last_S
    }
    
    if(print){
      setTxtProgressBar(pb, iter)
    }
  }
  
  save_res$execution_time = Sys.time() - start_time
  
  #output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = "empty",
  #                   all_graphs = niter-nburn, all_weights = all_weights, last_graph = last_G,
  #                   last_K = last_K, last_Theta = last_Theta )
  #
  
  if(print){close(pb)}
  return(save_res)
}

# NEW GIBBS OUTER ---------------------------------------------------------

Gibbs_sampler_update_h = function(
    set_UpdateParamsGSL_list,
    niter, 
    initialization_values_h, 
    alpha_target,
    alpha_add,
    adaptation_step,
    seed,
    update_sigma_prior,
    update_theta_prior, 
    update_weights,
    update_partition,
    update_graph,
    perform_shuffle,
    # new parameters for informed prior
    update_gamma,          # T or F
    rho_0,                 # initial expert partition
    eta,                    # prior probability for gamma
    compute_partition_update_info, # df for each iter for debugging,
    sample_eta,            # T or F
    # algorithm for graph
    algorithm_graph
){
  set.seed(seed)
  # Create a list for chains to save the values of each iteration
  chains <- list(
    Beta = vector("list", length = niter),
    mu = vector("list", length = niter),
    tau_eps = vector("list", length = niter),
    K = vector("list", length = niter),
    G = vector("list", length = niter),
    z = vector("list", length = niter),
    rho = vector("list", length = niter),
    time = vector("list", length = niter), 
    sigma = vector("list", length = niter),
    theta = vector("list", length = niter),
    weights_a = vector("list", length = niter),
    weights_d = vector("list", length = niter),
    H = vector("list", length = niter),
    gamma = vector("list", length = niter),
    partition_update_info = vector("list", length = niter),
    eta = vector("list", length = niter),
    how_many_accepted = vector("list", length = niter),
    # for graph
    bd_graph_start = vector("list", length = niter), 
    total_weights = vector("list", length = niter) 
  ) 
  
  # Initialization of the chains
  chains$Beta[[1]] <- initialization_values_h$Beta        
  chains$mu[[1]] <- initialization_values_h$mu
  chains$tau_eps[[1]] <- initialization_values_h$tau_eps 
  chains$K[[1]] <- initialization_values_h$K
  chains$G[[1]] <- initialization_values_h$G
  chains$z[[1]] <- initialization_values_h$z 
  chains$rho[[1]] <- initialization_values_h$rho 
  chains$time <- 0
  chains$sigma <- initialization_values_h$sigma 
  chains$theta <- initialization_values_h$theta
  chains$weights_a[[1]] <- initialization_values_h$weights_a
  chains$weights_d[[1]] <- initialization_values_h$weights_d
  
  chains$H[[1]] <- sum(initialization_values_h$gamma)  
  chains$gamma[[1]] <- initialization_values_h$gamma
  
  chains$eta[[1]] <- eta
  
  chains$bd_graph_start[[1]] <- initialization_values$graph_start 
  chains$total_weights[[1]] <- initialization_values$total_weights 
  
  # initialization of parameters for set_options
  weights_a <- initialization_values_h$weights_a
  weights_d <- initialization_values_h$weights_d
  total_weights <- initialization_values_h$total_weights       
  total_K <- initialization_values_h$total_K
  total_graphs <- initialization_values_h$total_graphs
  graph_start <- initialization_values_h$graph_start
  
  # gamma <- initialization_values_h$gamma not this, but like rho use chain$gamma[s-1]
  
  seeds = sample(1:99999, size = (niter + 10))
  
  pb = txtProgressBar(min = 2, max = niter, initial = 2, style = 3)  #initialize progress bar
  
  # print('Start Iters -----------------------------------------------------')
  for(s in 2:niter) {
    
    fit = UpdateParamsGSL(
      chains$Beta[[s-1]],
      chains$mu[[s-1]],
      chains$tau_eps[[s-1]], 
      chains$K[[s-1]],
      set_UpdateParamsGSL_list$tbase_base,
      set_UpdateParamsGSL_list$tbase_data,
      set_UpdateParamsGSL_list$Sdata,
      set_UpdateParamsGSL_list$a_tau_eps, 
      set_UpdateParamsGSL_list$b_tau_eps, 
      set_UpdateParamsGSL_list$sigma_mu,
      set_UpdateParamsGSL_list$r,
      set_UpdateParamsGSL_list$Update_Beta,
      set_UpdateParamsGSL_list$Update_Mu,
      set_UpdateParamsGSL_list$Update_Tau,
      seeds[s]
    )
    
    # Save Beta
    chains$Beta[[s]] <- fit$Beta
    
    # Save mu
    chains$mu[[s]] <- fit$mu
    
    # Save tau
    chains$tau_eps[[s]] <- fit$tau_eps
    
    # Set options for a single iteration of the Gibbs_sampler
    #print('Defining Options --------------------------------------------')
    options = set_options_h(
      sigma_prior_0 = chains$sigma[[s-1]],
      sigma_prior_parameters = list("a"=initialization_values_h$a_sigma,"b"=initialization_values_h$b_sigma,
                                    "c"=initialization_values_h$c_sigma,"d"=initialization_values_h$d_sigma),
      theta_prior_0 = chains$theta[[s-1]],
      theta_prior_parameters=list("c"=initialization_values_h$c_theta,"d"=initialization_values_h$d_theta),
      rho=chains$rho[[s-1]],    
      weights_a0 = weights_a,
      weights_d0 = weights_d,
      total_weights0 = total_weights,
      total_K0 = total_K,
      total_graphs0 = total_graphs,
      graph = graph_start,
      alpha_target = alpha_target,
      beta_mu = initialization_values_h$graph_density, 
      beta_sig2 = initialization_values_h$beta_sig2, 
      d = initialization_values_h$d, 
      alpha_add = alpha_add, 
      adaptation_step = adaptation_step,
      update_sigma_prior = update_sigma_prior,
      update_theta_prior = update_theta_prior,
      update_weights = update_weights,
      update_partition = update_partition,
      update_graph = update_graph,
      perform_shuffle = perform_shuffle,
      # informed partition
      gamma = chains$gamma[[s-1]],
      update_gamma = update_gamma,
      compute_partition_update_info = compute_partition_update_info,
      eta = chains$eta[[s-1]],
      sample_eta = sample_eta
    )
    
    # PRINT TO CANCEL ---------------------------------------------------------
    #print('Before options and after res done done')    
    #print('Gamma old in options that must be updated: '); print(chains$gamma[[s-1]])    
    
    # Run an iteration of the Gibbs Sampler
    res <- Gibbs_sampler_h(
      data = t(fit$Beta - fit$mu),
      niter = 1, nburn = 0, thin = 1,
      options = options,
      seed = seed,
      print = FALSE,
      s = s,         # Added this for weight adaptivity
      rho_0 = rho_0, # Added this for initial partition
      algorithm = algorithm_graph
    )
    
    z = do.call(rbind, lapply(res$rho, rho_to_z))
    
    # Save rho
    chains$rho[[s]] <- res$rho[[1]]
    
    # Save gamma and H and partition_update_info: informed partition
    chains$gamma[[s]] <- res$gamma[[1]]
    chains$H[[s]] <- sum(res$gamma[[1]])
    chains$partition_update_info[[s]] <- res$partition_update_info
    chains$how_many_accepted[[s]] <- res$accepted
    chains$eta[[s]] <- res$eta[[1]]
    # Save K
    chains$K[[s]] <- res$K [[1]]
    
    # Save G
    chains$G[[s]] <- res$G [[1]]    
    
    # Save z
    chains$z[[s]] <- z      
    
    # Save times for each K
    chains$time[[s]] <- res$execution_time
    
    # Save sigma and theta
    chains$sigma[[s]] <- res$sigma[[1]]
    chains$theta[[s]] <- res$theta[[1]]
    
    # Save weights
    chains$weights_a[[s]] <- res$weights_a[[1]]
    chains$weights_d[[s]] <- res$weights_d[[1]]
    
    # save graphs and respective weight
    chains$bd_graph_start[[s]] = res$bdgraph_start 
    chains$total_weights[[s]] = res$total_weights 
    
    # Update quantities for the next iteration
    weights_a <- res$weights_a[[1]]
    weights_d <- res$weights_d[[1]]
    total_weights <- res$total_weights
    total_K <- res$total_K[[1]]
    total_graphs <- res$total_graphs[[1]]
    graph_start = res$bdgraph_start
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  return(chains)
}


