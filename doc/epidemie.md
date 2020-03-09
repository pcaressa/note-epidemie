# Introduzione ai modelli matematici di diffusione epidemiologica
### Paolo Caressa, PhD
#### Ver. 1.0 marzo 2020

## 1. Premessa

Lo scopo di questi appunti è di fornire una rapida introduzione ai modelli deterministici della diffusione epidemica: si tratta delle nozioni più elementari, che presuppongono soltanto una conoscenza di base del calcolo differenziale, diciamo derivate in una variabile e opzionalmente integrali. Tuttavia spero che queste note diano una seppur fugace impressione dell'importanza dei metodi elementari dei sistemi dinamici, che veramente richiedono quattro derivate per essere espressi, e stimolino interesse scientifico nell'attuale congiuntura epidemica di cui tanto si (s)parla.

**Disclaimer:** Questi miei appunti sono stati scritti frettolosamente e potrebbero contenere errori e inesattezze, che prego di segnalare o correggere "forkando" il repository [https://www.github.com/pcaressa/note-epidemie](https://www.github.com/pcaressa/note-epidemie) sui quali li metto liberamente a disposizione. I programmi di simulazione sono stati scritti (in Python) non in modo ottimale ma per essere copiaincollati e modificati a piacimento, quindi risulteranno ai programmatori ineleganti e artigianali. Naturalmente non usateli per nessuna applicazione reale!!!

## 2. Il modello più semplice: nessuna cura...

Consideriamo il diffondersi di un virus in una popolazione: quel che ci interessa è la diffusione dell'epidemia, e per farlo suddividiamo la popolazione in "compartimenti".

Nel modello più semplice abbiamo due compartimenti:

- $S$ (che sta per *Suscettibili*) cioè il gruppo di persone che non hanno la malattia ma possono contrarla.
- $I$ (che sta per *Infetti*) cioè il gruppo di persone che hanno contratto la malattia.

Quel che vogliamo capire è il passaggio $S\to I$, che comporta capire quanti suscettibili divengono infetti, e a che velocità. In questo modello supponiamo quindi che un suscettibile possa divenire infetto e, una volta contagiato, permanga nel compartimento $I$.

Supporremo che il numero totale di persone nella popolazione sia $N$ e rimanga costante, mentre ovviamente supponiamo che $S$ e $I$ varino nel tempo, dunque ad ogni istante $t$ abbiamo che $S(t)+I(t)=N$. Quel che vogliamo determinare è l'incremento $\Delta I(t) = I(t+\Delta t) - I(t)$ di persone che contraggono l'infezione fra l'istante $t$ e l'istante $t+\Delta t$. L'ipotesi che facciamo è che questo incremento soddisfi la relazione di proporzionalità

$$
    \Delta I(t) = \beta \frac{I(t)S(t)}{N}\,\Delta t
$$

dove $\beta>0$ è una costante che denota di quanto crescono gli infetti nell'intervallo di tempo $\Delta t$, precisamente è il *tasso di contatti* nella popolazione. Supponendo che le nostre variabili siano continue (anche se in realtà sono discrete in quanto rappresentano la numerosità di una popolazione), siamo quindi tentati di dedurne che

$$
    \frac{\Delta I(t)}{\Delta t} = \beta \frac{I(t)S(t)}{N}
$$

e quindi, passando al limite per $\Delta t\to 0$, trovare la relazione differenziale

$$
    \frac{dI}{dt} = \beta \frac{IS}{N}
$$

Stante l'ipotesi che $N$ sia costante nel tempo, da cui $S(t)=N-I(t)$, possiamo anche scrivere l'equazione per la dinamica di $S$:

$$
    \frac{dS}{dt} = - \frac{dI}{dt} = -\frac\beta N IS
$$

A ben vedere queste due equazioni sono equivalenti all'unica equazione

$$
    I'(t) = \beta \frac{I(N-I)}{N}
$$

(Spesso si scrive $I'(t)$ per indicare la derivata $\frac{dI}{dt}$ di una funzione $I(t)$.)

Notiamo immediatamente che questa equazione esprime il fatto che la derivata di $I(t)$ è positiva (poiché lo è la costante $\beta$ e poiché lo sono gli altri termini che compaiono a secondo membro dell'equazione differenziale), e quindi ci aspettiamo che gli infettati crescano di continuo, fino a saturare l'intera popolazione, e svuotare il compartimento $S$ (a quel punto $I(N-I)=0$ e quindi avremmo da un certo istante in poi $I'=0$, cioè $I=N$ in modo costante).

Verifichiamo analiticamente e numericamente questa triste deduzione.


### Integrazione esatta del modello **SI**

Poiché la soluzione di questa equazione differenziale si può scrivere in modo esatto, facciamolo (non capita quasi mai!): questa parte è puramente matematica ma comunque esorto tutti a provare a leggerla in quanto ci consentirà di guardare "in faccia" la soluzione del nostro modello, cioè la funzione $I(t)$ che indica l'andamento degli infettati.

L'equazione che abbiamo ricavato, $I'=\beta I(N-I)/N$, ha la caratteristica che le variabili che vi figurano sono separabili, per cui possiamo portare tutti termini che contengono la $I$ da una parte lasciando quelli che contengono la $t$ dall'altra, integrando poi ambo i membri come (scrivo $I_0=I(t_0)$ e $I_1=I(t_1)$)

$$
    \int_{I_0}^{I_1} \frac{N}{I(N-I)}dI = \beta \int_{t_0}^{t_1} dt = \beta(t_1 - t_0)
$$

L'integrazione del primo membro di questa equazione è facile (uso la sostituzione $J=N-I$ da cui $dJ=-dI$):

\begin{align*}
    \int_{I_0}^{I_1} \frac{N}{I(N-I)}dI
        &= \int_{I_0}^{I_1} \frac{dI}{I} + \int_{I_0}^{I_1} \frac{dI}{N-I}
            = \int_{I_0}^{I_1} \frac{dI}{I} - \int_{N-I_0}^{N-I_1} \frac{dJ}{J}   \\
        &= \log\frac{I_1}{I_0} - \log\frac{N-I_1}{N-I_0}
            = \log\frac{I_1(N-I_0)}{I_0(N-I_1)}
\end{align*}

Pertanto abbiamo

$$
    e^{\beta(t_1-t_0)} = \frac{I_1(N-I_0)}{I_0(N-I_1)}
$$

da cui segue $I_0(N-I_1)e^{\beta(t_1-t_0)} = I_1(N-I_0)$, quindi $I_0Ne^{\beta(t_1-t_0)}=I_1(N-I_0+I_0e^{\beta(t_1-t_0)})$, il che ci consente di scrivere esplicitamente $I_1$:

$$
    I_1 = \frac{I_0Ne^{\beta(t_1-t_0)}}{N-I_0+I_0e^{\beta(t_1-t_0)}}
        = \frac{I_0N}{S_0e^{-\beta(t_1-t_0)}+I_0}
$$

rammentando che $S=I-N$.

Abbiamo quindi trovato una espressione esplicita per l'andamento degli infetti per ogni $t$,

$$
    I(t) = \frac{I_0N}{(N-I_0)e^{-\beta(t-t_0)}+I_0}
$$

che dipende da una serie di costanti: $I_0$ è il valore iniziale di $I$ mentre $\beta$ è una costante che dipende dalla malattia.

Proviamo a graficare questa funzione per alcuni valori di $\beta$, supponendo che la condizione iniziale, su una popolazione di $N=1000$ sia di $I_0=10$ infetti. Ecco un programmino Python che svolge questo compito.


```python
# Importiamo un po' di librerie che useremo qui e in seguito
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Dati del modello
N = 1000
t0 = 0
t1 = 10
t = np.linspace(t0, t1, 100)   # questo genera 100 valori fra t0 e t1

# Condizione iniziale = numero di infetti all'istante iniziale
I0 = 10

def grafico(handle, beta):
    y = [I0*N/((N-I0)*np.exp(-beta*(tt-t0))+I0) for tt in t]
    handle.plot(t, y, 'b', label=r'$\beta$=' + str(beta))
    handle.grid()
    handle.legend(loc='best')

# Imposta le aree per i 4 grafici e chiama la funzione che li disegna
# con vari valori dei parametri
fig, axes = plt.subplots(2, 2, sharex='all', sharey='all')

grafico(axes[0,0], beta = 0.5)
grafico(axes[0,1], beta = 1)
grafico(axes[1,0], beta = 2)
grafico(axes[1,1], beta = 10)

```


![png](output_3_0.png)


Come si vede, in ogni caso il numero degli infetti tende a crescere fino a saturare tutta la popolazione: il parametro $\beta$ agisce come un tasso di crescita, nel senso che maggiore è il suo valore e più in fretta la popolazione degli infetti tende a saturare l'intera popolazione.

Potete divertirvi a variare la condizione iniziale e $\beta$, ma alla fine il risultato sarà sempre quello.

### Integrazione numerica del modello **SI**

Il calcolo che abbiamo fatto nell'integrazione esatta è semplice e divertente: purtroppo so che qualcuno potrebbe opinare su questi giudizi, e anche in vista di esempi più complicati, mi sento quindi obbligato a fornire anche una integrazione dell'equazione differenziale con metodi numerici, cioè non basandoci sulle nostre magie algebriche simboliche ma lasciando fare al computer quello per cui è fatto, cioè i conti approssimati della soluzione del modello.

Si tratta di usare il solutore di equazioni differenziali ordinarie di Python e dargli in pasto direttamente la funzione che figura a secondo membro nell'equazione differenziale: ricordiamo che il modello **SI** è

$$
    I'(t) = \beta \frac{I(N-I)}{N}
$$

Quindi è del tipo $y' = f(t,y)$ (nel nostro caso $y=I$). Ecco un programmino che usa questo metodo, il quale non ci fornisce esplicitamente la funzione $I$ che cerchiamo, ma ce la calcola nei punti desiderati, che in fondo è tutto quel che ci serve.



```python
# Importiamo un po' di librerie che useremo qui e in seguito
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

N = 1000
t0 = 0
t1 = 10
t = np.linspace(t0, t1, 100)

# Condizione iniziale = numero di infetti all'istante iniziale
I0 = 10

def f(y,t, beta, N):
    return beta * y *(N - y)/N

def grafico(handle, beta):
    y = odeint(f, I0, t, args=(beta,N))    # risolve l'equazione differenziale
    handle.plot(t, y, 'b', label=r'$\beta$=' + str(beta))
    handle.grid()
    handle.legend(loc='best')

# Imposta le aree per i 4 grafici e chiama la funzione che li disegna
# con vari valori dei parametri
fig, axes = plt.subplots(2, 2, sharex='all', sharey='all')
grafico(axes[0,0], beta = 0.5)
grafico(axes[0,1], beta = 1)
grafico(axes[1,0], beta = 2)
grafico(axes[1,1], beta = 10)

```


![png](output_6_0.png)


Non è sorprendente che sia venuto lo stesso risultato precedente: l'integrazione numerica non è esatta ma nella pratica del calcolo scientifico fornisce gli stessi risultati di quella esatta (nella rappresentazione dei numeri all'interno dei computer non c'è nulla di realmente esatto...).

## Il modello **SIS**

I risultati del modello **SI** non sembrano molto incoraggianti, ma è chiaro che il modello è troppo semplificato rispetto alle epidemie reali, dato che la dinamica del modello prevede solo la possibilità di passare dal compartimento $S$ al compartimento $I$.

Nella realtà i malati possono anche guarire (per fortuna!), col che tornerebbero nel compartimento $S$: dunque un modello un pochino più sofisticato è dato dallo schema

$$
    S \to I \to S
$$

Questo passo ulteriore ignora alcune caratteristiche delle malattie infettive, cioè che ci si può immunizzare da esse o anche morire di esse, vedremo nella prossima sezione come includere anche queste caratteristiche nei nostri modelli, per ora continuiamo a fare un passo per volta.

Naturalmente in questo modello **SIS** le equazioni si complicano, per tenere conto anche della transizione $I\to S$, da infetto a suscettibile: in particolare, se nel modello **SI** l'equazione che descrive il tasso con cui cresce il comparimento degli infetti è

$$
    I'(t) = \beta \frac{I(N-I)}{N}
$$

nel modello **SIS** questo tasso di crescita viene frenato da una certa quota di infetti che tornano suscettibili, secondo un tasso di decrescita proporzionale a $I$ mercé un nuovo parametro $\gamma$:

$$
    I'(t) = \beta \frac{I(N-I)}{N} - \gamma I
$$

Questo parametro $\gamma$ si interpreta immaginando che $1/\gamma$ sia il tempo medio di durata dell'infezione in un individuo.

Dal punto di vista del compartimento dei suscettibili, l'equazione precedente diviene

$$
    S'(t) = -\beta\frac{S(N-S)}{N} +\gamma(N-S)
$$

E infatti ci torna: stiamo infatti dicendo che, ad ogni istante $t$, $S$ si decrementa per la percentuale di soggetti che passano dal compartimento $S$ al compartimento $I$, ma si incrementa per la percentuale dei soggetti che da $I$ tornano in $S$, e lo fanno con un tasso di incremento dato da $\gamma$.

Notiamo inoltre che, sommando queste equazioni, troviamo $I'+S' = 0$, cioè che $I+S$ è costante, cosa che supponevamo nel modello **SI**.

L'equazione per la dinamica della $I$ è ancora a variabili separabili e se qualcuno vuole divertirsi a dedurla come abbiamo fatto in precedenza gli dò il suggerimento di distinguere il caso $\beta=\gamma$ dal caso $\beta\neq\gamma$, poiché danno luogo ad equazioni diverse.

Qui cedo alla pigrizia e faccio fare il conto al computer: nei seguenti programmini Python mostro le curve relative alla soluzione $I(t)$ del modello **SIS** per diversi valori e combinazioni di $\beta$ e $\gamma$, pensando sempre a una popolazione di $N=1000$ persone con un numero iniziale di infetti $I_0=100$ (naturalmente cambiando i dati del programma si possono sperimentare le curve per altri valori, il motivo principale per cui queste note sono redatte in un foglio Jupyter è proprio dare a chi le legge la possibilità di giocare con i modelli).

Nei grafici indichiamo i valori del rapporto

$$
    \sigma = \frac{\beta}{\gamma}
$$



```python
# Importiamo un po' di librerie che useremo qui e in seguito
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

N = 1000
t0 = 0
t1 = 10
t = np.linspace(t0, t1, 100)

# Condizione iniziale = numero di infetti all'istante iniziale
I0 = 100

def f(y,t, beta, gamma, N):
    return beta * y *(N - y)/N - gamma * y

def grafico(handle, beta, gamma):
    y = odeint(f, I0, t, args=(beta,gamma,N))
    handle.plot(t, y[:,0], 'b', label="I(t)")
    handle.grid()
    handle.set_title(r'$\sigma$=' + str(beta/gamma))
    handle.legend(loc='best')

fig, axes = plt.subplots(2, 2, sharex='all', sharey='all')
grafico(axes[0,0], beta = 1, gamma = 0.1)
grafico(axes[0,1], beta = 1, gamma = 0.5)
grafico(axes[1,0], beta = 1, gamma = 1)
grafico(axes[1,1], beta = 1, gamma = 2)

```


![png](output_9_0.png)


Come si vede il numero degli infetti non tende necessariamente al numero totale della popolazione: in particolare, se $\sigma\leq 1$ il numero degli infetti decresce sempre di più, e in un tempo sufficientemente lungo si azzera.

Invece, se $\sigma>1$ il numero degli infetti cresce con un andamento simile a quello del modello SIS ma con l'importante differenza di non tendere a saturare tutta la popolazione, ma di assestarsi su un numero massimo di infetti che è percentualmente dato da $I_E = (\sigma-1)/\sigma$, e che si chiama *equilibrio endemico*: una volta raggiunto l'equilibrio endemico, il compartimento degli infetti si stabilizza e l'epidemia diviene una *endemia*, durevole nel tempo.

Per esempio, se $\sigma=2$ allora l'equilibrio endemico è raggiunto per un numero di infetti pari a $N\times I_E=N/2$. E infatti si vede dal grafico corrispondente che la curva raggiunge questo punto.

Per avere una idea della *asintoticità* di questo valore limite, ripetiamo il calcolo allungando smisuratamente i tempi:


```python
t0 = 0
t1 = 100
t = np.linspace(t0, t1, 1000)
fig, axes = plt.subplots(1, 1)
grafico(axes, beta = 1, gamma = 0.5)
```


![png](output_11_0.png)


Vediamo bene come l'equilibrio endemico costituisca una soglia che non è possibile superare. Per questo motivo la proprietà che stabilisce il comportamento asintotico del modello sulla base della soglia si suole chiamare **teorema della soglia critica**; nel caso del modello **SIS** è quello che abbiamo enunciato più sopra discutendo i casi $\sigma\leq1$ e $\sigma>1$.

PS: possiamo facilmente modificare il programmino per mostrare una animazione al variare di $\beta$: il risultato viene salvato in una gif animata...


```python
# Importiamo un po' di librerie che useremo qui e in seguito
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
#plt.style.use('seaborn-pastel')
import numpy as np
from scipy.integrate import odeint

N = 1000
t0 = 0
t1 = 10
t = np.linspace(t0, t1, 100)

# Condizione iniziale = numero di infetti all'istante iniziale
I0 = 100

def f(y,t, beta, gamma, N):
    return beta * y *(N - y)/N - gamma * y

fig = plt.figure()
ax = plt.axes(xlim=(t0, t1), ylim=(0, N))
line, = ax.plot([], [], lw=3)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    global beta
    y = odeint(f, I0, t, args=(beta,gamma,N))
    line.set_data(t, y)
    beta += 0.1
    return line,

beta = 0.1
gamma = 5
anim = FuncAnimation(fig, animate, init_func=init, frames=100, interval=1, blit=True)
anim.save('SIS.gif', writer='imagemagick')
anim

```




    <matplotlib.animation.FuncAnimation at 0x7f095558f208>




![png](output_14_1.png)


## Il modello **SIR**

Il modello **SIS** ha un che di ragionevole, ma anche di rassegnatamente drammatico: l'epidemia diviene una *endemia* (cioè il suo radicamento temporale si fa paragonabile alla vita degli individui: i modelli per le endemie utilizzano solitamente anche dei parametri che tengono conto di nascite e morti nella popolazione durante il tempo di diffusione della malattia).

Tuttavia è ragionevole supporre che il numero degli infetti diminuisca non solo perché alcuni di loro tornano nei suscettibili, ma anche per altri motivi: in generale questo accade per tutte le epidemie. Gli infetti possono guarire e sviluppare una immunità, dunque non essere più suscettibili di contrarre la malattia, o morire, o essere isolati.

Viene quindi spontaneo modificare il nostro modello introducendo un nuovo compartimento, $R$ (che sta per *rimossi*) e contiene gli infetti che non sono più tali ma nemmeno più suscettibili. Lo schema di questo modello **SIR** contempla che si possa passare nel compartimento $R$ solo se si proviene dal compartimento $I$:

$$
    S \to I \to R
$$

(Questa tipologia di modelli fu proposta da Kermack e McKendrick in una serie di classici lavori fra il 1927 e il 1933.)

In sostanza si tratta di aggiungere una equazione per la dinamica del passaggio $I\to R$, il che conduce ad ampliare le equazioni del modello **SI** come segue:

$$
\begin{cases}
\displaystyle
    S' = -\beta\frac{IS}{N}   \\
\displaystyle
    I' = \beta \frac{IS}{N} -\gamma I  \\
\displaystyle
    R' = \gamma I
\end{cases}
$$

Infatti stiamo dicendo che la quota $\gamma I$ che viene sottratta dal compartimento $I$ a ogni $\Delta t$ passa nel compartimento $R$, che quindi si incrementa di $\gamma I\Delta t$.

Notiamo che stavolta $N=S+I+R$ e che, sommando le tre equazioni differenziali, troviamo $S'+I'+R'=0$, dunque di nuovo $N$ costante. In particolare possiamo sempre determinare $R(t)$ una volta che siano note $S(t)$ e $I(t)$, quindi possiamo limitarci a studiare le prime due equazioni per avere la dinamica dei tre compartimenti.

Quel che viene fuori è che anche in questo caso quel che governa è la costante $\sigma = \beta/\gamma$, e che il teorema della soglia critica per il modello **SIR** afferma quanto segue.

- Se $\sigma\leq1$ allora $I(t)$ decresce indefinitamente.
- Se $\sigma>1$ e $S_0 \leq N/\sigma$ allora $I(t)$ decresce indefinitamente.
- Se $\sigma>1$ e $S_0 > N/\sigma$ allora $I$ cresce fino a un certo istante $T$ e poi decresce indefinitamente.

Come si vede è interessante l'ultima: c'è un tempo $T$ corrispondente al **picco dell'epidemia**, il punto di massimo della funzione $I(t)$...

Una spiegazione intuitiva è la seguente: il massimo di una funzione si ha in corrispondenza di dove si annulla la sua derivata, e dato che la derivata di $I(t)$ è il secondo membro dell'equazione differenziale $I'(t)=\beta SI/N -\gamma I$, troviamo che il punto di massimo di $I$ si ha per l'istante $T$ tale che $S(T) = N\gamma/\beta = N/\sigma$.

Ora, se $\sigma\leq 1$ $I(t)$ decresce sicuramente, altrimenti cresce: ma se $S_0\leq N/\sigma$, dato che inizialmente $S$ decresce (in quanto cresce $I$), $I$ non potrà assumere il suo massimo, e quindi decrescerà indefinitamente. Se invece all'istante $t_0$ si ha che $S_0>N/\sigma$, al crescere di $I(t)$ si avrà un decrescere di $S(t)$ che, supponendola continua, prima o poi dovrà toccare il valore soglia di $N/\sigma$ e quindi far giungere $I(t)$ al suo massimo: da lì in poi, $I(t)$ non può che decrescere.

Una illustrazione di questo risultato è data dal seguente programmino, che risolve il sistema di equazioni differenziali in $S$ e $I$, e deduce i valori di $R$ come $N-S-I$.



```python
# Importiamo un po' di librerie che useremo qui e in seguito
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

N = 1000
t0 = 0
t1 = 10
t = np.linspace(t0, t1, 100)

# Condizione iniziale = numero di infetti all'istante iniziale
S0 = 900
I0 = 100

def f(y,t, beta, gamma, N):
    return [-beta * y[0] *y[1]/N, beta * y[0] *y[1]/N - gamma * y[1]] 

def grafico(handle, beta, gamma):
    y = odeint(f, [S0,I0], t, args=(beta,gamma,N))
    handle.plot(t, y[:,1], 'b', label="I(t)")
    handle.plot(t, N - y[:,0] - y[:,1], 'r', label="R(t)")
    handle.grid()
    handle.legend(loc='best')
    handle.set_title(r'$\sigma$=' + str(beta/gamma))

fig, axes = plt.subplots(2, 2, sharex='all', sharey='all')
grafico(axes[0,0], beta = 1, gamma = 0.1)
grafico(axes[0,1], beta = 1, gamma = 0.5)
grafico(axes[1,0], beta = 1, gamma = 2)
grafico(axes[1,1], beta = 1, gamma = 10)

```


![png](output_16_0.png)


Si vede chiaramente nei grafici più in alto che se $\sigma>1$ la funzione $I(t)$ cresce fino a un certo tempo $T$ e poi decresce. Facendo un conto non troppo complicato, si trova che

$$
    I(T) = I_0 - \frac{N}{\sigma}\left(1+\log\frac{\sigma S_0}{N}\right)
$$

(Si supponga che $I=F(S(t))$ per una certa funzione $F$: derivando questa definizione e usando le due equazioni differenziali come definizioni di $S'(t)$ e $I'(t)$, si giunge all'equazione precedente, tenendo conto che $S(T)=N/\sigma$.)

Ultima osservazione: la $\sigma$ nel modello **SIR** corrisponde al famoso coefficiente $R_0$ di cui ogni tanto si leggono stime e al quale si attribuisce molta importanza, giustificata dal teorema della soglia critica.

## Ulteriori modelli

Il modello **SIR** sembra più realistico dei precedenti in quanto contempla la possibilità che dalla malattia si guarisca, si muoia o ci si immunizzi. E non propone invariabilmente uno scenario apocalittico di endemia.

Prima di applicare questo modello a dei calcoli concreti, vale la pena di aggiungere che fin qui abbiamo scalfito la punta dell'iceberg: infatti esistono due modi di generalizzare questo modello.

1. Aggiungere ulteriori compartimenti, come per esempio gli *esposti*, cioè chi ha contratto la malattia ma non è ancora contagioso, se non dopo un periodo di latenza. Quindi il modello diviene $S\to E\to I\to R$. Un altro compartimento si ha introducendo il compartimento $M$ dei neonati che hanno una immunità temporanea, che può durare qualche mese, indotta dagli anticorpi della madre: in questo caso il modello diviene $M\to S\to E\to I\to R$. Naturalmente a ogni compartimento si aggiunge una equazione differenziale e un parametro, quindi le cose tendono a farsi "multidimensionali" e complicate, tant'é che non credo sia noto nessun teorema di soglia critica per il modello **MSEIR**.

2. Estendere il modello epidemico a un modello endemico che tenga conto di nascite e morti durante l'evoluzione del contagio: in questo caso non si aggiungono altre equazioni ma altri parametri. Per esempio **MSEIR** è inerentemente endemico, considerando anche i neonati.

Per chi fosse interessato ad approfondire queste tematiche, come primo riferimento in italiano per la parte teorica di questo affascinante argomento suggerisco l'ultimo capitolo del manuale di Mascia-Montefusco *Un invito alla biomatematica*, LaDotta, 2015.

Un testo semplice e accessibile che offre una introduzione alla matematica delle malattie infettive ma anche una panoramica dei metodi di stima parametrica necessari per "incrociarli" con i dati reali è il recente libro di Li *An Introduction to Mathematical Modeling of Infectious Diseases* Springer, 2018.

A livello più avanzato trovo ancora illuminante e assolutamente consigliabile la trattazione matematica di Hethcote *The Mathematics of Infectious Diseases* SIAM Review, 42 (2000), 599-653.


## Facciamo lavorare i modelli

Per quelli come me che trovano piacere estetico nella manipolazione di formule e nella deduzione dei teoremi, quanto ho scritto in queste note basta per solleticare qualche neurone. Ma, specie a fronte dell'attuale emergenza da Covid 19, non posso chiudere questi appunti senza tentare una simulazione usando dati reali.

Tengo a sottolineare come questo sia un semplice esercizio per mostrare l'uso di un modello, non una proposta di predizione. Userò il modello **SIR**, il che vuol dire che dobbiamo prendere i dati di una epidemia e in particolare, per ciascun istante, il totale dei suscettibili, degli infetti e dei rimossi.

E ora vengono i problemi: qualcuno ha efficacemente detto che *viviamo in una stitichezza di informazioni ma sommersi da una diarrea di dati*. Girando in rete si trovano moltissimi dati su questa influenza, e tanti grafici basati su modelli con le quattro operazioni (numero di morti nell'unità di tempo, rapporto fra eposti e infetti, etc.). Ma ho fatto fatica a trovare i dati utili per il modello.

### Un dataset e come utilizzarlo

Ai miei fini servono compartimenti misurati, quindi userò come popolazione totale quella dei sottoposti a un controllo di presenza del virus (anche se questo dato è drogato dal fatto che normalmente si considerano i sintomatici per fare i controlli). Gli infetti sono quelli risultati positivi. I rimossi sono sia i morti sia chi è guarito, supponendo (cosa non completamente vera a quanto pare) che un rimosso non possa tornare suscettibile.

Nei dataset normalmente si elencano per paese i casi positivi, senza indicare rispetto a che campione della popolazione, un dato ai nostri fini inutile.

Comunque su GitHub ci sono alcuni repository con dati che fanno al caso nostro: in particolare ne ho trovato uno della Protezione Civile Italiana, che consiglio a tutti quelli che vogliano attingere ai dati sul Covid19 in Italia, e precisamente al repository [https://github.com/pcm-dpc/COVID-19](https://github.com/pcm-dpc/COVID-19) gestito da Umberto Rosini. In particolare troviamo una serie storica in formato JSON, che ho copiaincollato alla data in cui sto scrivendo queste righe nello script qui di seguito, che la memorizza in una variabile DATA:


```python
DATA = [
    {
        "data": "2020-02-24 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 101,
        "terapia_intensiva": 26,
        "totale_ospedalizzati": 127,
        "isolamento_domiciliare": 94,
        "totale_attualmente_positivi": 221,
        "nuovi_attualmente_positivi": 221,
        "dimessi_guariti": 1,
        "deceduti": 7,
        "totale_casi": 229,
        "tamponi": 4324
    },
    {
        "data": "2020-02-25 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 114,
        "terapia_intensiva": 35,
        "totale_ospedalizzati": 150,
        "isolamento_domiciliare": 162,
        "totale_attualmente_positivi": 311,
        "nuovi_attualmente_positivi": 90,
        "dimessi_guariti": 1,
        "deceduti": 10,
        "totale_casi": 322,
        "tamponi": 8623
    },
    {
        "data": "2020-02-26 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 128,
        "terapia_intensiva": 36,
        "totale_ospedalizzati": 164,
        "isolamento_domiciliare": 221,
        "totale_attualmente_positivi": 385,
        "nuovi_attualmente_positivi": 74,
        "dimessi_guariti": 3,
        "deceduti": 12,
        "totale_casi": 400,
        "tamponi": 9587
    },
    {
        "data": "2020-02-27 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 248,
        "terapia_intensiva": 56,
        "totale_ospedalizzati": 304,
        "isolamento_domiciliare": 284,
        "totale_attualmente_positivi": 588,
        "nuovi_attualmente_positivi": 203,
        "dimessi_guariti": 45,
        "deceduti": 17,
        "totale_casi": 650,
        "tamponi": 12014
    },
    {
        "data": "2020-02-28 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 345,
        "terapia_intensiva": 64,
        "totale_ospedalizzati": 409,
        "isolamento_domiciliare": 412,
        "totale_attualmente_positivi": 821,
        "nuovi_attualmente_positivi": 233,
        "dimessi_guariti": 46,
        "deceduti": 21,
        "totale_casi": 888,
        "tamponi": 15695
    },
    {
        "data": "2020-02-29 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 401,
        "terapia_intensiva": 105,
        "totale_ospedalizzati": 506,
        "isolamento_domiciliare": 543,
        "totale_attualmente_positivi": 1049,
        "nuovi_attualmente_positivi": 228,
        "dimessi_guariti": 50,
        "deceduti": 29,
        "totale_casi": 1128,
        "tamponi": 18661
    },
    {
        "data": "2020-03-01 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 639,
        "terapia_intensiva": 140,
        "totale_ospedalizzati": 779,
        "isolamento_domiciliare": 798,
        "totale_attualmente_positivi": 1577,
        "nuovi_attualmente_positivi": 528,
        "dimessi_guariti": 83,
        "deceduti": 34,
        "totale_casi": 1694,
        "tamponi": 21127
    },
    {
        "data": "2020-03-02 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 742,
        "terapia_intensiva": 166,
        "totale_ospedalizzati": 908,
        "isolamento_domiciliare": 927,
        "totale_attualmente_positivi": 1835,
        "nuovi_attualmente_positivi": 258,
        "dimessi_guariti": 149,
        "deceduti": 52,
        "totale_casi": 2036,
        "tamponi": 23345
    },
    {
        "data": "2020-03-03 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 1034,
        "terapia_intensiva": 229,
        "totale_ospedalizzati": 1263,
        "isolamento_domiciliare": 1000,
        "totale_attualmente_positivi": 2263,
        "nuovi_attualmente_positivi": 428,
        "dimessi_guariti": 160,
        "deceduti": 79,
        "totale_casi": 2502,
        "tamponi": 25856
    },
    {
        "data": "2020-03-04 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 1346,
        "terapia_intensiva": 295,
        "totale_ospedalizzati": 1641,
        "isolamento_domiciliare": 1065,
        "totale_attualmente_positivi": 2706,
        "nuovi_attualmente_positivi": 443,
        "dimessi_guariti": 276,
        "deceduti": 107,
        "totale_casi": 3089,
        "tamponi": 29837
    },
    {
        "data": "2020-03-05 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 1790,
        "terapia_intensiva": 351,
        "totale_ospedalizzati": 2141,
        "isolamento_domiciliare": 1155,
        "totale_attualmente_positivi": 3296,
        "nuovi_attualmente_positivi": 590,
        "dimessi_guariti": 414,
        "deceduti": 148,
        "totale_casi": 3858,
        "tamponi": 32362
    },
    {
        "data": "2020-03-06 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 2394,
        "terapia_intensiva": 462,
        "totale_ospedalizzati": 2856,
        "isolamento_domiciliare": 1060,
        "totale_attualmente_positivi": 3916,
        "nuovi_attualmente_positivi": 620,
        "dimessi_guariti": 523,
        "deceduti": 197,
        "totale_casi": 4636,
        "tamponi": 36359
    },
    {
        "data": "2020-03-07 18:00:00",
        "stato": "ITA",
        "ricoverati_con_sintomi": 2651,
        "terapia_intensiva": 567,
        "totale_ospedalizzati": 3218,
        "isolamento_domiciliare": 1843,
        "totale_attualmente_positivi": 5061,
        "nuovi_attualmente_positivi": 1145,
        "dimessi_guariti": 589,
        "deceduti": 233,
        "totale_casi": 5883,
        "tamponi": 42062
    }
]
```

Ai miei fini i dati che mi interessano sono i valori di $I$ e $R$, che posso estrarre come

1. $I(t)$ è dato dal valore del campo "totale_attualmente positivi";
2. $R(t)$ è dato dalla somma dei valori dei campi "dimessi_guariti" e "deceduti"
3. $N$ è dato dal valore del campo "tamponi".

Rispetto ai nostri modelli notiamo come il valore $N$ sia tutto meno che costante: un modo di ovviare a questo inconveniente è prendere il valore più recente (e quindi più alto) e usarlo come $N$ anche per tutti gli istanti precedenti.

Ma in realtà vorremmo svincolarci da questo $N$, e per farlo possiamo modificare il nostro modello per trattare valori relativi e non assoluti.

## Una estensione del modello **SIR**

Le nostre equazioni differenziali, vale la pena ricordarle ancora una volta, sono

$$
\begin{cases}
\displaystyle
    S' = -\beta\frac{IS}{N}   \\
\displaystyle
    I' = \beta \frac{IS}{N} -\gamma I  \\
\displaystyle
    R' = \gamma I
\end{cases}
$$

Tuttavia nella nostra serie storica, notiamo come il valore $N$ della popolazione totale (che sto assumendo essere il numero di persone sottoposte a tampone) non sia costante: dobbiamo quindi trasformare $N$ in una variabile, e l'ipotesi più comune (sulla cui giustificazione non mi soffermerò qui) è che segua una legge esponenziale:

$$
    \frac{dN}{dt} = \alpha N
$$

dove $\alpha$ è il tasso di crescita della popolazione. Attenzione: *nel nostro caso il tasso di crescita della popolazione non è dovuto alla differenza fra nuove nascite e nuovi decessi nell'unità di tempo, ma nell'aumento della popolazione per cui è stato fatto il test di infezione del vaccino (il famoso tampone)*.

Una volta resa $N$ una variabile, dobbiamo estendere il modello per tenere conto della crescita della popolazione: per farlo modifichiamo le equazioni come segue

$$
\begin{cases}
\displaystyle
    S' = -\beta\frac{IS}{N} + \alpha N   \\
\displaystyle
    I' = \beta \frac{IS}{N} -\gamma I \\
\displaystyle
    R' = \gamma I
\end{cases}
$$

Infatti i nuovi "censiti" nella nostra popolazione confluiscono a priori nel compartimento dei suscettibili, non degli infetti o dei rimossi. Aggiungendo l'equazione differenziale per l'incremento della popolazione troviamo il sistema

$$
\begin{cases}
\displaystyle
    S' = -\beta\frac{IS}{N} + \alpha N   \\
\displaystyle
    I' = \beta \frac{IS}{N} -\gamma I \\
\displaystyle
    R' = \gamma I   \\
\displaystyle
    N' = \alpha N
\end{cases}
$$

Effettivamente ora $S'+I'+R'=N'$ come ci aspettiamo. Se $\alpha=0$ ritroviamo il modello precedente.


### Calibrare il modello

Torniamo al nostro dataset: abbiamo un insieme di 13 misurazioni che ci consentono di ottenere, per 13 giorni consecutivi, i valori $I$, $R$ e $N$, da cui otteniamo $S$ per differenza. Quel che ci preme è chiaramente determinare $I$ e $R$, quindi concentrarci sulla seconda e terza equazione.

Naturalmente, prima di applicarlo a una situazione reale, dobbiamo trovare dei valori specifici per i parametri $\alpha, \beta, \gamma$ dai quali dipende il modello stesso, e i cui valori abbiamo visto incidere in modo decisivo sulla forma delle soluzioni.

Cominciamo da $\alpha$: in questo caso abbiamo 13 rilevazioni per $N(t)$, con $t=0,1,2,...,12$, precisamente


```python
serie_N = [misura["tamponi"] for misura in DATA]
serie_N
```




    [4324,
     8623,
     9587,
     12014,
     15695,
     18661,
     21127,
     23345,
     25856,
     29837,
     32362,
     36359,
     42062]



Poiché ci aspettiamo che questi siano i valori di $N(0), N(1), ..., N(12)$ possiamo determinare l'$\alpha$ che renda la soluzione dell'equazione $N'=\alpha N$ nei punti dati il più vicina possibile ai valori dati dalla serie storica.

Questa operazione di il parametro di una equazione differenziale che renda la soluzione più coerente con un insieme di valori dati si chiama *stima parametrica*, e rientra nella smisurata tematica dei metodi di ottimizzazione: non aggiungerò altro a questo se non che la libreria SciPy di Python ce ne mette a disposizione quanti ne vogliamo.

Per non complicare troppo le cose (e lasciare il piacere di implementare i minimi quadrati a chi legge) userò il seguente criterio: consideriamo $\alpha$ una variabile e calcoliamo con un metodo numerico il minimo rispetto a questa variabile della funzione ottenuta prendendo la norma (cioè la radice delle somme dei quadrati delle differenze) fra i vettori dati dai punti della serie storica e dai valori desunti risolvendo l'equazione differenziale per $\alpha$. Questa funzione dipende da $\alpha$, e minimizzandola rispetto a questa variabile troviamo il valore di $\alpha$ per cui commettiamo (in media quadratica) l'errore minimo nell'usare il modello invece che i dati reali.

Come si vede dal seguente programmino, che usa la funzione fmin che minimizza una funzione data senza stare troppo a chiedersi come sia fatta, è più facile farlo che dirlo.


```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fmin

serie_N = np.array([misura["tamponi"] for misura in DATA])
N0 = serie_N[0]

t0 = 0
t1 = len(serie_N)
intervallo = np.linspace(t0, t1, len(serie_N))

def f(y, t, alpha):
    return alpha*y[0]

def modello(alpha):
    y = odeint(f, N0, intervallo, args=(alpha,))
    return np.linalg.norm(y[:,0] - serie_N)

# Trova il valore ottimale di alpha
alpha_ottimale = fmin(modello, 2.5)

print("alpha ottimale =", alpha_ottimale[0])

# Ora disegna la soluzione per l'alpha ottimale
y = odeint(f, N0, intervallo, args=(alpha_ottimale,))
plt.plot(np.linspace(t0, t1, 13), serie_N, 'ro', label="Dati reali")
plt.plot(intervallo, y[:,0], label="Modello")
plt.legend()
```

    Optimization terminated successfully.
             Current function value: 16548.147606
             Iterations: 26
             Function evaluations: 52
    alpha ottimale = 0.18399429321289062
    




    <matplotlib.legend.Legend at 0x23bc89f5508>




![png](output_27_2.png)


Sembra quindi che il tasso di crescita nell'effettuare i tamponi possa assumersi come $\alpha=0.184$ circa. Poiché la soluzione dell'equazione $N'=\alpha N$ è

$$
    N = N_0e^{\alpha (t-t_0)} = 4324 e^{0.184 t}
$$

le rimanenti equazioni dipendono solo da $S$ e $I$: in effetti dobbiamo stimare contemporaneamente i parametri $\beta$ e $\gamma$ usando le due equazioni per $S'$ e $I'$ simultaneamente, che è quel che facciamo nel prossimo programmino.


```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fmin

serie_I = np.array([misura["totale_attualmente_positivi"] for misura in DATA])
I0 = serie_I[0]

serie_R = np.array([misura["dimessi_guariti"] + misura["deceduti"] for misura in DATA])
R0 = serie_R[0]

serie_N = np.array([misura["tamponi"] for misura in DATA])
N0 = serie_N[0]

serie_S = serie_N - serie_I - serie_R
S0 = serie_S[0]

serie_SIR = np.array([[serie_S[i], serie_I[i], serie_R[i]] for i in range(len(DATA))])
SIR0 = serie_SIR[0]

t0 = 0
t1 = len(serie_N)
intervallo = np.linspace(t0, t1, len(serie_N))

alpha_ottimale = 0.184

def fun_N(t):
    return N0 * np.exp(alpha_ottimale*(t-t0))

def fun_SIR(y, t, beta, gamma):
    # y = [S(t), I(t), R(t)]
    return [-beta*y[0]*y[1]/fun_N(t)+alpha_ottimale*fun_N(t),
            beta*y[0]*y[1]/fun_N(t)-gamma*y[1],
            gamma*y[1] ]

def modello(betagamma):
    # Risolve il sistema di equazioni $S'=..., I'=..., R'=...$
    y = odeint(fun_SIR, SIR0, intervallo, args=(betagamma[0],betagamma[1]))
    return np.linalg.norm(y - serie_SIR)

# Trova il valore ottimale di beta e gamma
valori_ottimali = fmin(modello, (1, 1))

print("beta ottimale =", valori_ottimali[0])
print("gamma ottimale =", valori_ottimali[1])

# Ora disegna la soluzione per l'alpha ottimale
y = odeint(fun_SIR, SIR0, intervallo, args=(valori_ottimali[0],valori_ottimali[1]))
plt.plot(np.linspace(t0, t1, 13), serie_I, 'ro', label="Dati reali per I")
plt.plot(intervallo, y[:,1], label="Modello per I")
plt.legend()
```

    Optimization terminated successfully.
             Current function value: 15827.343280
             Iterations: 39
             Function evaluations: 74
    beta ottimale = 0.29471246460464284
    gamma ottimale = 0.015319160219554274
    




    <matplotlib.legend.Legend at 0x23bc8b32f88>




![png](output_29_2.png)


Come si vede il modello è ora calibrato abbastanza bene sui dati reali: questo non garantisce sulla sua predittività, anche perché $N$ varia non come una popolazione ma come un test somministrato a una popolazione.

Per curiosità tracciamo comunque le curve $S$, $I$ e $R$ su un intervallo temporale di 30 giorni.


```python
intervallo = np.linspace(t0, t0 + 30)
y = odeint(fun_SIR, SIR0, intervallo, args=(valori_ottimali[0],valori_ottimali[1]))
plt.plot(intervallo, y[:,0], "g", label="Modello per S")
plt.plot(intervallo, y[:,1], "b", label="Modello per I")
plt.plot(intervallo, y[:,2], "r", label="Modello per R")
plt.legend()
```




    <matplotlib.legend.Legend at 0x23bc8ba3ec8>




![png](output_31_1.png)


Dunque, il modello prevede una impennata di infezioni che crescono indefinitamente, questo perché supporre una crescita esponenziale in una popolazione per tempi medio-lunghi non è realistico, si tratta di una modellizzazione che  funziona solo nel breve periodo a meno di non bilanciare il tasso di crescita con un tasso di decrescita legato a decessi, etc.: un ulteriore modo di affinare il modello sarebbe sostituire $N$ con una funzione che si comporta inizialmente come l'esponenziale ma poi tende asintoticamente a un certo valore che è quello presunto dei tamponi che saranno fatti in totale.

Ma in ogni caso lo scopo di questa modellizzazione è stato vedere come il modello può spiegare l'andamento attuale: lanciarsi in analisi predittive richiede una continua ricalibrazione perché i $\beta$ e $\gamma$ (che determinano il comportamento asintotico del modello) in realtà cambiano nel tempo.


## Conclusioni

Non è stato lo scopo di questa breve nota quello di proporre un modello predittivo per le epidemie, né tantomeno predire alcunché riguardo l'attuale Covid19.

Piuttosto lo scopo è stato mostrare:

1. come i modelli matematici sono di complessità incrementale, e crescente e come operino delle ipotesi comunque semplicistiche sulla natura delle cose;
2. che il modello senza i dati o i dati senza modello non dicono nulla;
3. che i modelli hanno dei limiti di applicabilità e che nell'applicare il modello stesso si operano delle scelte, dettate dalla semplicità, dagli strumenti a disposizione, etc.;
4. che alla fine oggi è molto facile programmare un modello, rappresentare i risultati in forma grafica e giocare con i parametri dai quali il modello dipende.

In particolare, molti lettori di queste note avranno immediatamente visto un migliore uso del modello che ho proposto, o una sua modifica per renderlo più efficiente, o reperire più dati per calibrarlo meglio, etc. etc.

Ma se il messaggio che non basta esporre una tabella o un grafico per stabilire alcunché su un fenomeno complesso come une epidemia è passato, lo scopo di questa breve nota (oltre al divertimento nello scriverla) sarà stato raggiunto.
