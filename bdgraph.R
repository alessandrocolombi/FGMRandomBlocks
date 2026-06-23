# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C)  Laura Codazzi, Alessandro Colombi, Matteo Gianella                           |
#                                                                                                 |
#     Part of this code is adapted from BDgraph package (C Reza Mohammadi)                        |
#                                                                                                 |
#     FGM is free software: you can redistribute it and/or modify it under                        |
#     the terms of the GNU General Public License as published by the Free                        |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
#                                                                                                 |
#     Maintainer: Laura Codazzi (laura.codazzi@tuhh.de),                                          |
#                 Alessandro Colombi (a.colombi10@campus.unimib.it),                              |
#                 Matteo Gianella (matteo.gianella@polimi.it)                                     |
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


# Old function ------------------------------------------------------------

#' BDgraph function for sampling from GGM
#' @param g.prior May take a single value, a full pxp matrix or a vector of length 2. The first two cases refer to a Bernoulli prior where the parameter is equal for all the possible links (first case)
#' or it is link specific (second case). The third case refers to the multiplicity correction prior (Bernoulli-Beta) and it takes the two Beta hyperparameters.
#' @export
bdgraph = function( data, rho, n = NULL, method = "ggm", algorithm = "bdmcmc", iter = 5000,
                    burnin = iter / 2, not.cont = NULL, g.prior = 0.5, df.prior = 3,
                    CCG_D = NULL, g.start = "empty", jump = NULL, save = FALSE, print = 1000,
                    cores = NULL, threshold = 1e-8, beta_params )
{
    if( df.prior < 3  ) stop( " 'prior.df' must be >= 3" )
    if( iter < burnin ) stop( " Number of iteration must be more than number of burn-in" )
    burnin <- floor( burnin )
    if( print > iter ) print = iter
    
    #Converting the representation of the current partition from rho to z
    z = rho_to_z(rho)
    groups_cardinality=rho
    n_groups=length(rho)

    cores = get_cores( cores = cores )

    list_S_n_p = get_S_n_p( data = data, method = method, n = n, not.cont = not.cont )
    S      = list_S_n_p $ S
    n      = list_S_n_p $ n
    p      = list_S_n_p $ p
    method = list_S_n_p $ method
    colnames_data = list_S_n_p $ colnames_data

    if( method != "ggm" ) stop('FGM package supports only ggm method.')

    D = diag(p)
    b      = df.prior
    b_star = b + n
    Ds     = D + S
    Ts     = chol( solve( Ds ) )
    Ti     = chol( solve( D ) )   # only for double Metropolic-Hastings algorithms

    g_prior = get_g_prior( g.prior = g.prior, p = p )
    FGM_gprior_length = dim(g_prior)[1]*dim(g_prior)[2]

    G       = get_g_start( g.start = g.start, g_prior = g_prior, p = p )
    K       = get_K_start( G = G, g.start = g.start, Ts = Ts, b_star = b_star, threshold = threshold )

    Theta=get_S_from_G_rho(G, rho)

    if( save == TRUE ){
        qp1           = ( p * ( p - 1 ) / 2 ) + 1
        string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
        sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100"
        graph_weights = c( rep ( 0, iter - burnin ) )         # waiting time for every state
        all_graphs    = c( rep ( 0, iter - burnin ) )         # vector of numbers like "10100"
        all_weights   = c( rep ( 1, iter - burnin ) )         # waiting time for every state
        size_sample_g = 0
    }else{
        p_links = matrix( 0, p, p )
    }

    #if( ( save == TRUE ) && ( p > 50 & iter > 20000 ) )
    #{
    #    cat( "  WARNING: Memory needs to run this function is around " )
    #    print( ( iter - burnin ) * utils::object.size( string_g ), units = "auto" )
    #}

    K_hat      = matrix( 0, p, p )
    last_graph = K_hat
    last_K     = K_hat

    if( ( is.null( jump ) ) && ( p > 10 & iter > ( 5000 / p ) ) ){
        jump = floor( p / 10 )
    }

    if( is.null( jump ) ) jump = 1

    if( ( p < 10 ) && ( jump > 1 ) )      cat( " WARNING: the value of jump should be 1 " )
    if( jump > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: the value of jump should be smaller " )
    # mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
    # cat( mes, "\r" )

    # - -  main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - - |
    if( save == TRUE ){
        
      if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) ){
            result = .C( "ggm_bdmcmc_map",as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior),  as.integer(FGM_gprior_length),
                       as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                       all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat),
                       sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                       as.integer(b), as.integer(b_star), as.double(Ds), as.integer(print), Theta = as.integer(Theta), z = as.integer(z), as.integer(n_groups),
                       groups_cardinality = as.integer(groups_cardinality), 
                       beta_params = as.double(beta_params), PACKAGE = "FGM" )
      }else{
            stop('FGM can handle only the jump == 1 case. Code should be cleaned up.')
      }
    }else{
      stop('FGM package can not handle the save==FALSE option in the bdgraph function')
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    K_hat      = matrix( result $ K_hat, p, p, dimnames = list( colnames_data, colnames_data ) )
    last_graph = matrix( result $ G    , p, p, dimnames = list( colnames_data, colnames_data ) )
    last_K     = matrix( result $ K    , p, p )
    last_Theta = matrix( result $ Theta, n_groups, n_groups)

    if( save == TRUE ){
        
      if( algorithm == "rjmcmc" ) K_hat = K_hat / ( iter - burnin )
      
      size_sample_g = result $ size_sample_g
      sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
      graph_weights = result $ graph_weights[ 1 : size_sample_g ]
      all_graphs    = result $ all_graphs + 1
      all_weights   = result $ all_weights
      
      if( ( algorithm != "rjmcmc" ) & ( jump != 1 ) )
      {
          all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
          all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ]
      }

      output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat,
                     all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph,
                     last_K = last_K, last_Theta = last_Theta )
    }else{
      stop('FGM package can not handle the save==FALSE option in the bdgraph function')
    }

    class( output ) = "bdgraph"
    return( output )
  }




# New function ------------------------------------------------------------

#' post_graph_sampling function for sampling from GGM
#' @param g.prior May take a single value, a full pxp matrix or a vector of length 2. The first two cases refer to a Bernoulli prior where the parameter is equal for all the possible links (first case)
#' or it is link specific (second case). The third case refers to the multiplicity correction prior (Bernoulli-Beta) and it takes the two Beta hyperparameters.
#' @export

post_graph_sampling = function( data, rho, n = NULL, method = "ggm", algorithm = "rjmcmc", iter = 5000,
                    burnin = iter / 2, not.cont = NULL, g.prior = 0.5, df.prior = 3,
                    CCG_D = NULL, g.start = "empty", jump = NULL, save = TRUE, print = 1000,
                    cores = NULL, threshold = 1e-8, beta_params ){
    if( df.prior < 3  ) stop( " 'prior.df' must be >= 3" )
    if( iter < burnin ) stop( " Number of iteration must be more than number of burn-in" )
    burnin <- floor( burnin )
    if(is.null(print) || !is.finite(print) || print <= 0)
      print = iter
    if(print > iter)
      print = iter
    print = as.integer(max(1L, print))
    
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data must be a matrix or dataframe" )
    if( is.data.frame( data ) ) data <- data.matrix( data )
    if( any( is.na( data ) ) ) stop( "the algorithm does not deal with missing values" )	
    

    #Converting the representation of the current partition from rho to z
    z = rho_to_z(rho)
    groups_cardinality=rho
    n_groups=length(rho)

    cores = get_cores( cores = cores )

    list_S_n_p = get_S_n_p( data = data, method = method, n = n, not.cont = not.cont )
    S      = list_S_n_p $ S
    n      = list_S_n_p $ n
    p      = list_S_n_p $ p
    method = list_S_n_p $ method
    colnames_data = list_S_n_p $ colnames_data

    # p = ncol(data)
    # n = nrow(data)
    # S <- t(data) %*% data
    # colnames_data = NULL
    
    # any(S!=list_S_n_p $ S)
    # class(S)
    # class(list_S_n_p $ S)
    # 
    # n==list_S_n_p$n
    # class(n)
    # class(list_S_n_p $ n)
    # 
    # p==list_S_n_p$p
    # class(p)
    # class(list_S_n_p $ p)
    
    if( method != "ggm" ) stop('FGM package supports only ggm method.')

    D      = diag(p)
    b      = df.prior
    b_star = b + n
    Ds     = D + S
    Ts     = chol( solve( Ds ) )
    Ti     = chol( solve( D ) )   # only for double Metropolic-Hastings algorithms

    g_prior = get_g_prior( g.prior = g.prior, p = p )
    FGM_gprior_length = dim(g_prior)[1]*dim(g_prior)[2]

    G = get_g_start( g.start = g.start, g_prior = g_prior, p = p )
    K = get_K_start( G = G, g.start = g.start, Ts = Ts, b_star = b_star, threshold = threshold )

    Theta = get_S_from_G_rho(G, rho)

    if( save == TRUE ){
        qp1           = ( p * ( p - 1 ) / 2 ) + 1
        string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
        sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100"
        graph_weights = c( rep ( 0, iter - burnin ) )         # waiting time for every state
        all_graphs    = c( rep ( 0, iter - burnin ) )         # vector of numbers like "10100"
        all_weights   = c( rep ( 1, iter - burnin ) )         # waiting time for every state
        size_sample_g = 0
    }else{
      stop('FGM package cannot handle the save==FALSE option in post_graph_sampling() function')
    }

    K_hat      = matrix( 0, p, p )
    last_graph = K_hat
    last_K     = K_hat

    if( ( is.null( jump ) ) && ( p > 10 & iter > ( 5000 / p ) ) ){
        jump = floor( p / 10 )
    }

    if( is.null( jump ) ) jump = 1

    if( ( p < 10 ) && ( jump > 1 ) )      cat( " WARNING: the value of jump should be 1 " )
    if( jump > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: the value of jump should be smaller " )
    
    # - -  main BDMCMC & RJMCMC algorithms implemented in C++ - - - - - - - - - - - - - |
    if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) ){
      
      result = .C( "ggm_bdmcmc_map", 
                   as.integer(iter), 
                   as.integer(burnin), 
                   G = as.integer(G), 
                   as.double(g_prior),  
                   as.integer(FGM_gprior_length),
                   as.double(Ts), 
                   K = as.double(K), 
                   as.integer(p), 
                   as.double(threshold),
                   all_graphs = as.integer(all_graphs), 
                   all_weights = as.double(all_weights), 
                   K_hat = as.double(K_hat),
                   sample_graphs = as.character(sample_graphs), 
                   graph_weights = as.double(graph_weights), 
                   size_sample_g = as.integer(size_sample_g),
                   as.integer(b), 
                   as.integer(b_star), 
                   as.double(Ds), 
                   as.integer(print), 
                   Theta = as.integer(Theta), 
                   z = as.integer(z), 
                   as.integer(n_groups),
                   groups_cardinality = as.integer(groups_cardinality), 
                   beta_params = as.double(beta_params), 
                   PACKAGE = "FGM" )
      
    }else if(( method == "ggm" ) && ( algorithm == "rjmcmc" )){
      
      result = .C( "ggm_rjmcmc_map",
                   as.integer(iter), 
                   as.integer(burnin), 
                   G = as.integer(G), 
                   as.double(g_prior),  
                   as.double(Ts), 
                   K = as.double(K), 
                   as.integer(p), 
                   as.double(threshold),
                   all_graphs = as.integer(all_graphs), 
                   all_weights = as.double(all_weights), 
                   K_hat = as.double(K_hat),
                   sample_graphs = as.character(sample_graphs), 
                   graph_weights = as.double(graph_weights), 
                   size_sample_g = as.integer(size_sample_g),
                   as.integer(b), 
                   as.integer(b_star), 
                   as.double(Ds), 
                   as.integer(print), 
                   Theta = as.integer(Theta), 
                   z = as.integer(z), 
                   as.integer(n_groups),
                   groups_cardinality = as.integer(groups_cardinality), 
                   beta_params = as.double(beta_params), 
                   PACKAGE = "FGM" )
      
    }else if(( method == "ggm" ) && ( algorithm == "rjmcmc.mpl" )){
      
      result = .C( "ggm_rjmcmc_mpl_map", 
                   as.integer(iter), 
                   as.integer(burnin), 
                   G = as.integer(G), 
                   as.double(g_prior), 
                   as.double(S), 
                   as.integer(n), 
                   as.integer(p), 
                   all_graphs = as.integer(all_graphs), 
                   all_weights = as.double(all_weights), 
                   sample_graphs = as.character(sample_graphs), 
                   graph_weights = as.double(graph_weights), 
                   size_sample_g = as.integer(size_sample_g),
                   as.integer(print), 
                   Theta = as.integer(Theta), 
                   z = as.integer(z), 
                   as.integer(n_groups),
                   groups_cardinality = as.integer(groups_cardinality),
                   beta_params = as.double(beta_params), 
                   PACKAGE = "FGM" )
      
    }else{
      
    stop('FGM package supports only algortihm == "rjmcmc" or "rjmcmc.mpl" or "bdmcmc"')
    
    }
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    
    size_sample_g = result $ size_sample_g
    sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
    graph_weights = result $ graph_weights[ 1 : size_sample_g ]
    all_graphs    = result $ all_graphs + 1
    all_weights   = result $ all_weights
    last_graph = matrix( result $ G, p, p )
    last_Theta = matrix( result $ Theta, n_groups, n_groups)
    
    if( algorithm == "rjmcmc" | algorithm == "bdmcmc" ) {
      K_hat = matrix( result $ K_hat, p, p )
      K_hat = K_hat / ( iter - burnin )
    }
    
    if(algorithm == 'rjmcmc.mpl'){
      
      output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, 
                     all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph,
                     last_Theta = last_Theta) 
      # .mpl does not give precision matrix, hence:
      #output$last_K = BDgraph::rgwish(adj = output$last_graph, b = b_star, D = Ds)
      output$last_K = get_K_start( G = output$last_graph, g.start = output$last_graph, 
                                   Ts = Ts, b_star = b_star, threshold = threshold )
     
      }else{
       
       last_K = matrix( result $ K, p, p )
       
       if( ( algorithm != "rjmcmc" ) & ( jump != 1 ) ){
         all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
         all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ]
        }
       
       output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat,
                     all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph,
                     last_K = last_K, last_Theta = last_Theta )
      }

    class( output ) = "bdgraph"
    return( output )
}














