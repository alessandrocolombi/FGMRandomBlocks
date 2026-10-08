# PriorK_simulation

## Simulazioni: soltanto la tabella finale

Dalla root del progetto sulla VM:

```bash
./run_job.sh PriorK_simulation/PriorK_simulation.R choose_wd=2 eta_usr=0.5
./run_job.sh PriorK_simulation/PriorK_sim_ab_theta.R choose_wd=2 eta_usr=0.5
```

`choose_wd=1` seleziona `C:/Users/colom/FGMRandomBlocks`;
`choose_wd=2` seleziona `/home/colombi/FGMRandomBlocks`.
Senza choose_wd resta la directory corrente. Le opzioni possono essere date
in qualsiasi ordine. eta_usr deve essere finito e in [0,1]. Il default e' 1;
per Ctrl+Invio si modifica la variabile nella configurazione iniziale.

Gli eta sono costruiti come richiesto:

```r
eta <- numeric(40)
eta[c(4, 6, 9, 13, 18, 22, 28, 33)] <- eta_usr
eta[40] <- 1
```

L'ultima assegnazione mantiene il vincolo del modello: il nodo 40 deve sempre
chiudere il dominio. eta_usr riguarda quindi effettivamente gli otto punti
esperti interni. Con eta_usr=0 c'e' un unico sottogruppo; con eta_usr=1 le
nove dimensioni sono (4,2,3,4,5,4,6,5,7). Ai valori intermedi i sottogruppi
vengono ricampionati a ogni replica.

- `PriorK_simulation.R`: griglia theta/sigma fissati, scrive soltanto
  `results/PriorK_mean_sd_table.csv` (righe theta, colonne sigma).
- `PriorK_sim_ab_theta.R`: sigma~Beta(1,1), theta+sigma~Gamma(a_theta,b_theta),
  con b_theta RATE; griglie a_theta/b_theta. Scrive soltanto
  `results/PriorK_ab_theta_mean_sd_table.csv`.

Entrambi stampano solo la tabella finale con celle media (sd) a due decimali.
Le combinazioni di parametri fuori dominio sono NA. Gli errori effettivi
(numerici, di compilazione o di input) restano segnalati; non vengono nascosti.
La sd e' quella di K, non l'errore Monte Carlo della media.
Non si salvano estrazioni, configurazioni o altri riepiloghi.
I vecchi output gia' presenti non vengono cancellati. I CSV mantengono i nomi
precedenti: una nuova esecuzione, anche con eta_usr diverso, li sovrascrive.
Il launcher Bash continua a mostrare PID e percorso del log per gestire il job.

Per source(), caricare lo script e poi chiamare esplicitamente la funzione:

```r
source("PriorK_simulation/PriorK_sim_ab_theta.R")
eta_usr <- 0.5
eta <- make_prior_eta(eta_usr)
simulate_prior_K_ab_theta(a_theta_grid, b_theta_grid, eta, B, seed, output_file)
```

## Grafico riproducibile di K

`Plot_PriorK_ab_theta.R` si esegue dall'alto con Ctrl+Invio, dalla cartella
progetto oppure PriorK_simulation. Modificare a_theta, b_theta, a_sigma, b_sigma, eta_usr, B e seed
nella configurazione. Usa lo stesso campionatore della tabella gerarchica.
Le barre sono le probabilita' Monte Carlo di K; non una KDE continua.

Le figure sono salvate in `results/figures/` con nomi come:

```
PriorK_a_theta_0.88_b_theta_0.94_eta_usr_0.5_a_sigma_1_b_sigma_1.pdf
PriorK_a_theta_0.88_b_theta_0.94_eta_usr_0.5_a_sigma_1_b_sigma_1.png
```

Nessun CSV o campione viene salvato dal grafico. `save_figures=FALSE` mostra
soltanto il grafico nel pannello Plots. Ripetere gli stessi parametri sovrascrive
le relative figure.

`draw_K_hyperprior(..., a_sigma=1, b_sigma=1)` permette di personalizzare la
Beta di sigma. Entrambi i parametri devono essere positivi e finiti; i default
preservano Beta(1,1) e il comportamento della tabella gerarchica.

## Coefficienti e verifiche

`PriorK_helpers.R` condivide parsing, eta e campionamento senza side effect.
`PriorKh_test.R` fornisce PMF e media verificate; `log_C_PY.cpp` calcola in log
la ricorrenza di Lijoi, Mena e Pruenster (2007):

    C_mathcal(n,k) = (n-1-sigma*k) C_mathcal(n-1,k) + sigma C_mathcal(n-1,k-1).

La PMF usa C_mathcal/sigma^k senza fattori di segno esterni. La variante scalata
comprende sigma=0 esattamente. I valori grezzi vengono verificati per massa,
intervallo e media (11) prima del campionamento.
Il vecchio compute_logC.cpp, lasciato intatto, usa n-1+sigma*k quando chiamato
con scale=-sigma e quindi non fornisce i coefficienti richiesti per sigma>0.

Richiesti R, Rcpp e un compilatore C++ (Rtools su Windows).
I test e gli script diagnostici mantengono i propri output espliciti; la
restrizione alla sola tabella si applica ai due script di simulazione.
