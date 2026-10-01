# Prior di K: implementazione corretta e verificata

## Risultato

La PMF e la media (11) sono coerenti. Il problema era il coefficiente restituito
chiamando il vecchio C++ con scale=-sigma, non il prefattore della PMF.
Il file originale compute_logC.cpp e il bridge esistente sono conservati.
La nuova implementazione e' in log_C_PY.cpp.

Per p=40, theta=1, sigma=0.5, senza normalizzazione aggiuntiva:

- massa: 1.000000000000014 (percorso con coefficienti scalati);
- media numerica: 12.40631636137310;
- formula (11): 12.40631636137297;
- tutte le 40 probabilita' sono tra 0 e 1.

La valutazione diretta log_C - k*log(sigma), nel controllo interattivo, differisce
solo per arrotondamento (errore della media circa 1.4e-13).

## Convenzione ed errore precedente

Lijoi, Mena e Pruenster (2007), p.774 e appendice p.782, usano

    mathcal_C(n,k;sigma) = 1/k! sum_{j=0}^k (-1)^j choose(k,j) (-j*sigma)_n,

con fattoriale crescente. Per 0<sigma<1 questi coefficienti sono positivi
per 1<=k<=n. La PMF e'

    P(K_n=k) = prod_{i=1}^{k-1}(theta+i*sigma)
               / (sigma^k * (theta+1)_{n-1}) * mathcal_C(n,k;sigma).

L'appendice specifica mathcal_C(n,k;sigma)=(-1)^(n-k)*C(n,k;sigma),
con C nella convenzione di Charalambides. Questa relazione di segno e' corretta.
Tuttavia il valore restituito da compute_logC(n,-sigma,0) NON e' il logaritmo
positivo di mathcal_C(n,k;sigma). La ricorrenza effettiva del vecchio codice,
ponendo s=sigma>0 e A=exp(output), e'

    A(n,k) = (n-1+sigma*k)*A(n-1,k) + sigma*A(n-1,k-1).

Essa corrisponde invece a (-1)^k*mathcal_C(n,k;-sigma). In particolare,
A(2,1)=sigma*(1+sigma), mentre mathcal_C(2,1;sigma)=sigma*(1-sigma).
Quindi non basta correggere il segno esterno: cambiano i valori assoluti.
La precedente risposta che applicava (-1)^n assumeva la corrispondenza dichiarata
tra output e C; il riferimento consente ora di identificarne precisamente
l'incompatibilita' con il codice. Con la convenzione mathcal_C del riferimento
non va aggiunto alcun (-1)^k o (-1)^n alla PMF.

Esempio theta=1, sigma=0.5, n=2:

    corretto: P(K=1)=0.25, P(K=2)=0.75, somma=1;
    vecchio:  valori=0.75 e 0.75, somma=1.5.

Nel Draft, se C indica quello di Charalambides nella convenzione esplicitata
nell'appendice del riferimento, il fattore di conversione e' (-1)^(n-k),
non (-1)^k. Usare direttamente mathcal_C elimina l'ambiguita'.
La formula (11) della media non necessita di correzioni.

## Nuova implementazione e controlli

La ricorrenza corretta, derivabile dalla definizione dell'appendice, e'

    mathcal_C(n,k) = (n-1-sigma*k)*mathcal_C(n-1,k)
                     + sigma*mathcal_C(n-1,k-1),
    mathcal_C(0,0)=1, mathcal_C(n,0)=0 per n>=1.

Il C++ usa log-sum-exp. Internamente calcola S=mathcal_C/sigma^k:

    S(n,k) = (n-1-sigma*k)*S(n-1,k) + S(n-1,k-1).

Questo evita cancellazioni quando sigma e' piccolo e comprende esattamente
sigma=0 (Stirling non segnati, limite Dirichlet). log_C_PY restituisce
log(mathcal_C) per k=0,...,n; log_C_scaled_PY restituisce log(S).
L'indice R corrispondente a k e' k+1.

I controlli sono eseguiti sui valori grezzi: nessuna divisione per la loro somma.
PriorKh_test.R confronta:

- massa, intervallo e media (11), a p=40 per nove coppie di parametri;
- coefficienti con la formula esplicita per sigma=1/2 dell'appendice,
  per tutti gli n=1,...,40 (errore massimo sui log circa 5.7e-14);
- PMF con una ricorrenza indipendente basata sulle probabilita' predittive PY;
- casi limite con sigma=0, sigma vicino a 0 o 1, theta negativo ma >-sigma.

## Esecuzione

Richiede R, Rcpp e un compilatore C++ (Rtools su Windows).
Impostare la working directory alla cartella del progetto o PriorK_simulation.
I percorsi sono semplici e gli script possono essere eseguiti con Ctrl+Invio.
Dalla cartella principale, in R:

    source("PriorK_simulation/PriorKh_test.R")
    run_prior_Kh_tests()
    source("PriorK_simulation/PriorKh_normalization_check.R")
    source("PriorK_simulation/PriorK_simulation.R")
    simulate_prior_K(theta_grid, sigma_grid, eta, B, seed)

Con Rscript ciascuno script esegue direttamente il proprio controllo/simulazione.
Il controllo interattivo confronta i nuovi coefficienti con il C++ originale.
Il file storico results/PriorKh_weights_diagnostic.csv conserva i vecchi pesi
errati; le PMF corrette sono in PriorKh_probabilities.csv e
PriorKh_normalization_values.csv.

## Simulazione di K

La sezione CONFIGURAZIONE UTENTE in PriorK_simulation.R contiene theta_grid,
sigma_grid, eta, B=1000 e il seme. p resta 40. Dominio: 0<=sigma<1,
theta>-sigma; quando sigma=0 occorre theta>0. Le combinazioni non ammissibili vengono saltate: statistiche e celle della
tabella sono NA, status=invalid_parameters e nessuna estrazione viene generata
(campi delle estrazioni NULL nel file RDS). Le altre coppie sono simulate normalmente.
I controlli numerici sulle PMF delle coppie valide rimangono attivi.

Gli eta richiesti sono 1 per j=4,6,9,13,18,22,28,33,40 e 0 altrove.
Quindi gamma=eta quasi certamente, H=9 e p_tilde=(4,2,3,4,5,4,6,5,7).
Per altre sequenze eta si campionano Bernoulli indipendenti; eta[40]=1 e'
obbligatorio. Gli eta non esperti devono essere zero per rispettare il Draft.
Condizionatamente a gamma si estraggono K_h indipendenti dalle rispettive PMF,
poi K=sum(K_h). Per gli eta attuali E[K]=sum_h E[K_h] e' deterministico.

Eseguite 1000 repliche per ciascuna delle nove combinazioni delle griglie
iniziali theta=(0,1,5), sigma=(0.2,0.5,0.8). I file in results includono:

- PriorK_grid_summary.csv: media, varianza campionaria, mediana, quantili
  0.025/0.975 (type=1), estremi, errore Monte Carlo della media e media
  delle aspettative condizionate;
- CSV per coppia con theta e sigma nel nome;
- RDS per coppia con configurazione, eta, griglie, seed, sessione, K, H,
  gamma, p_tilde, K_h e medie condizionate;
- PriorK_eta.csv con la sequenza usata.

Rieseguire nella stessa cartella sovrascrive i file con lo stesso nome.
La riproducibilita' richiede anche lo stesso ordine delle griglie.

## Avvio sulla VM

Dalla root del progetto sulla VM:

```bash
./run_job.sh PriorK_simulation/PriorK_simulation.R choose_wd=2
```

`choose_wd=1` seleziona `C:/Users/colom/FGMRandomBlocks`;
`choose_wd=2` seleziona `/home/colombi/FGMRandomBlocks`.
I percorsi sono modificabili in `wd_vec` all'inizio dello script.
Senza argomento viene mantenuta la working directory corrente.
Il launcher inoltra gli argomenti a Rscript e salva il log in `logs/`
relativamente alla directory da cui viene lanciato. Il percorso dello script
passato al launcher deve esistere prima della selezione della working directory.

L'output conclusivo e' `PriorK_mean_sd_table.csv`: theta sulle righe, sigma
sulle colonne, ciascuna cella nel formato `media (sd)` a due decimali.
La stessa tabella e' stampata a console (quindi nel log della VM).
La sd e' la deviazione standard campionaria di K, non l'errore della media.
Il riepilogo completo resta in `PriorK_grid_summary.csv`.
