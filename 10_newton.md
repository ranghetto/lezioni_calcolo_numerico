## 10. Metodo di Newton: convergenza locale, ordine di convergenza, test di arresto, esempi, altri metodi di linearizzazione

Partiamo dalla relazione chiave ottenuta nel capitolo precedente per l'errore del metodo di Newton:
$$
\begin{array}{c}
e_{n+1}=c_ne_n^2 \leq ce_n^2 \\
c_n = \dfrac{1}{2}\dfrac{|f''(z_n)|}{|f'(x_n)|},\quad c=\dfrac{1}{2}\dfrac{M_2}{m_1} \\
M_2 = \underset{x\in[c,d]}{\max} |f''(x)|, \; m_1 = \underset{x\in[c,d]}{\min} |f'(x)|\gt0
\end{array}
$$
assumendo che il metodo sia convergente, $f\in C^2[a,b]$ e che ${x_n}\subset[c,d]\subseteq[a,b]$ con $f'(x)\neq0\;\forall\;\in[c,d]$.

Da questa abbiamo ottenuto, fissato $\Theta\in(0,1)$ e preso $\bar{n}$ tale che $ce_n \leq \Theta\lt 1\;\forall\;n\ge\bar{n}$, $e_{\bar{n}+k}\le \dfrac{1}{c}\Theta^{2^k}, k\ge 0$ che ci ha fatto capire che per zeri semplici il metodo converge *piu' che esponenzialmente*.

A questo punto possiamo chiederci cosa succede se gia' con la scelta iniziale $x_0$ vale $ce_0\lt1$.
In questo caso avremmo la disuguaglianza
$$ce_n\le(ce_0)^2, n\ge0$$
infatti
$$\begin{array}{r}
ce_1 \le (ce_0)^2 \\
ce_2 \le (ce_1)^3 \le (ce_0)^4 \\
\vdots \\
ce_n\le(ce_{n-1})^2 \le (ce_0)^{2^n}
\end{array}$$
Questo ci fa intuire che se $ce_0 \lt1$ cioe' $e_0 = |x_0-\xi|\lt \dfrac{1}{c}$ cioe' se prendiamo $x_0$ in un intorno opportuna di $\xi$ zero di $f$ avremo la convergenza con le sole ipotesi che $f\in C^2$ e che $f'(x)\neq0$ in quell'intorno, perche' $ce_0\lt1 \implies (ce_0)^{2^n}\to0, n\to\infty \implies e_n\to0, n\to\infty$.

Si puo' infatti dimostrare il seguente teorema:

> **Teorema** (convergenza *locale* del metodo di Newton)
> Sia $\xi$ zero di $f$ ed $\exists\delta\gt0: f\in C^2(I_\delta)$ e $f'(x)\neq0\;\forall\;x\in I_\delta$, dove $I_\delta = [\xi-\delta, \xi+\delta]$; inoltre sia $x_0 \in (\xi-\gamma, \xi+\gamma)$ con $\gamma=\min\{\delta, \frac{1}{c}\}$ dove $c = \frac{1}{2}\dfrac{\underset{x\in I_\delta}{\max}|f''(x)|}{\underset{x\in I_\delta}{\min|f'(x)|}}$
> $$\implies\forall\;n\ge0,\;x_n\in I_\gamma, \; e_n\le\dfrac{1}{c}(ce_0)^{2^n}\to0, n\to\infty$$

Come abbiamo gia' osservato nella lezione precedente, il metodo di Newton pur essendo molto veolce sarebbe poco utile se funzionasse solo in ipotesi forti come ad esempio quelle del teorema di convergenza globale, in cui $e_0$ puo' essere grande ma le ipotesi su $f$ sono molto restrittive ($f''(x)$ di segno costante in $[a,b]$).

Invece, nel teorema di convergenza locale (che di nuovo, vale la pena di ribadirlo e' un set di condizioni sufficienti per la convergenza) basta che $f\in C^2$ e $f'\neq0$ in un intorno di $\xi$, con la seconda ipotesi che e' sempre soddisfatta se lo *zero e' semplice* perche' $f'(x)$ e' continua (permanenza del segno di $f'(\xi)$).

Conviene a questo punto formalizzareil concetto di *convergenza quadratica* del metodo di Newton, dando delle definizioni generali sull'ordine di convergenza di un metodo.

> **Definizione** (ordine di convergenza)
> Dato un metodo che produce una successione $\{x_n\}_{n\ge0}$ tale che $\lim_{n\to\infty} x_n=l$ con $l$ limite finito si dice che:
> 1. il metodo ha ordine di convergenza *almeno* $p\ge1$ se $\exists\;c\gt0$ (con $c \in (0,1)$ se $p=1$) tale che $e_{n+1}\ge ce_n^p \;\forall\;n$;
> 2. il metodo ha ordine di convergenza *esattamente* $p\ge1$ se $\exists\;L\gt0$ (con $L\in(0,1)$ se $p=1$) tale che $\lim_{n\to\infty}\dfrac{e_{n+1}}{e_n^p}=L$ (dove $L$ e' spesso chiamata *costante asintotica* del metodo) la convergenza e' detta **lineare** se $p=1$, **superlineare** se $p\gt1$ e in particolare **quadratica** se $p=2$, **cubica** se $p=3$, $\dots$

Facciamo alcune osservazioni: nel caso $p=1$ le condizioni $e_{n+1}\le ce_n, 0\lt c\lt1$ (convergenza almeno lineare) e $\lim_{n\to\infty}\dfrac{e_{n+1}}{e_n}=L$ con $L\in(0,1)$ sono condizioni sufficienti per la convergenza.
Infatti da $e_{n+1}\le ce_n$ si ricava $e_1\le ce_0, e_2\le ce_1\le c^2e_0, e_3\le ce_2 \le c^3e_0, \dots e_n\le c^ne_0$ con $c^n\to0,n\to\infty$ perche' $c\in(0,1)$ da cui $e_n\to 0, n\to\infty$.

Analogamente, se $\lim_{n\to\infty}\dfrac{e_{n+1}}{e_n}=L$ con $L\in(0,1)$, preso $\bar{\epsilon}\in(0,1-L)\;\exists\;\bar{n}$ tale che $0\le\dfrac{e_{n+1}}{e_n}\le L+\bar{\epsilon}=c\lt1 \;\forall\;n\ge\bar{n}$ cioe' $e_{n+1}\le ce_n \;\forall\;n\ge\bar{n}$, da cui ragionando come sopra $e_n\le c^{n-\bar{n}}e_\bar{n}\to0, n\to\infty$.
D'altra parte, se c'e' convergenza, allora necessariamente $L\le1$. Infatti, se $L\gt1$ allora l'errore divergerebbe, cioe' non ci sarebbe convergenza.

Per quanto riguarda il emtodo di Newton, dalla definizione di ordine di convergenza possiamo dire che l'ordine e' almeno $p=2$ per zeri semplici perche' sappiamo che in tal caso $e_{n+1}\le ce_n^2$; d'altra parte dall'analisi fatta sappiamo che
$$\dfrac{e_{n+1}}{e_n^2}= c_n = \dfrac{1}{2}\dfrac{|f''(z_n)|}{|f'(x_n)|}\to \dfrac{1}{2}\dfrac{|f''(\xi)|}{|f'(\xi)|}$$ per $n\to\infty$ (visto che $z_n\in\ \text{int}(\xi, x_n)$ e $f', f''$ sono continue quindi lo sono $|f'|,|f''|$ e si puo' portare il limite "dentro le funzioni").
Ma allora l'ordine e' esattamente $p=2$ se $f''(\xi)\neq0$, altrimenti $L=0$ e si puo' dimostrare che l'ordine e' almento $p=3$ se $f$ ammette derivata terza in $\xi$.

Vale la pena di notare che il metodo di bisezione, pur comportandosi "in media" come un metodo di ordine $p=1$ e costante asintotica $L=\frac{1}{2}$, non ha un ordine definito.
Infatti, in generale non si puo' dimostare che $\exists\;c:e_{n+1}\le ce_n$ con $0\lt c\lt1$ e tanto meno che $\lim_{n\to\infty}\dfrac{e_{n+1}}{e_n}=\dfrac{1}{2}$.

Facciamo infine un'ultima osservazione sull'ordine di convergenza del metodo di Newton: cosa succede quando $f'(\xi)=0$ cioe' quando lo zero non e' semplice.
In questo caso l'ordine di convergenza del metodo di Newton scende a $p=1$ con costante asintotica $1-\dfrac{1}{m}$ dove $m$ e' la molteplicita' dello $0$ cioe
$$f(\xi)=f'(\xi)=\cdots=f^{(m-1)}(\xi)=0 \; \text{e} \; f^{(m)}(\xi)\neq0$$
Ad esempio se $m=2$ cioe' se $f'(\xi)=0$ e $f''(\xi)\neq0$ si ha
$$\dfrac{e_{n+1}}{e_n}=\dfrac{e_{n+1}}{e_n^2}\cdot e_n=c_n\cdot e_n = \dfrac{1}{2}\dfrac{|f''(z_n)|}{|f'(x_n)|}\cdot e_n$$
Usando la formula di Taylor $f'(x_n)=\underset{=0}{f'(\xi)}+f''(\xi)(x_n-\xi)+o(e_n)$, da cui $\dfrac{|f'(x_n)|}{e_n}\to|f''(\xi)|, n\to\infty$ e quindi $\dfrac{e_{n+1}}{e_n}=e_nc_n\to\dfrac{1}{2},n\to\infty$.

Dopo questo "excursus" sull'importante nozione di ordine di convergenza, applicata al metodo di Newton, vediamo di spiegare come mai nel caso del calcolo di $\sqrt2$ il numero di cifre corrette essenzialmente raddoppia ad ogni iterazione.

Quando si parla di cifre corrette, che sono cifre di mantissa, si sta parlando sostanzialmente dell'errore relativo $r_n=\dfrac{e_n}{|\xi|}, \xi\neq0$.
Ora, da $e_{n+1}\le ce_n^2$ otteniamo
$$r_{n+1} = \dfrac{e_{n+1}}{|\xi|}\le c \dfrac{e_n^2}{|\xi|}=c|\xi|r_n^2$$
Se $c|\xi|\le1$ si ha $r_{n+1}\le r_n^2$, cioe' l'errore *relativo* al passo $n+1$ e' maggiorato dal quadrato dell'errore relativo al passo $n$: $r_1\le r_0^2, r_2\le r_1^2\le r_0^4,\dots, r_n\le r_{n-1}^2 \le r_0^{2^n}$.

Osserviamo che le richieste computazionali del metodo di Newton sono molto piu' forti di quelle del metodo di bisezione.
Se la bisezione infatti, nella versione base, ha bisogno solo di conoscere il segno di $f(x_n)$ (per cui basta come sappiamo un errore $\lt100\%$ sulle quantita') e col test del residuo pesato (che assume zero semplice e $f\in C^1$) basta saper stimare l'ordine di grandezza del residuo e del peso, ovvero l'ordine di grandezza di $f$ ed $f'$ almeno per $n$ abbastanza grande (infatti per la derivata possiamo accontentarci di uuna stima tramite rapporti incrementali), col metodo di Newton diventa essenziale poter calcolare ocn accuratezza sia $f(x_n)$ sia $f'(x_n)$ perche' gli errori su queste quantita tendono ad esseere dominanti nel processo iterativo di calcolo.

D'altra parte per il metodo di Newton non abbiamo una stima a priori facile da usare per arrestare le iterazioni.
Pensiamo infatti ad esempio al caso della convergenza locale e alla stima $e_n\le\dfrac{1}{c}(ce_0)^{2^n}$.
Questa ci dice qualitativamente che c'e' convergenza prendendo $x_0$ in un intorno di $\xi$, ma per trovare tale intorno e per usare la stima avremmo bisogno di conoscere $c$, che richiede di saper stimare il minimo di $|f'|$ ma anche il massimo di $|f''|$ (il che e' agevole quando di $f$ e' nota un'espressione analitica semplice ma non in generale).

Per fortuna, il metodo di Newton puo' sfruttare una semplice stima a posteriori, il cosiddetto **step** $|x_{n+1}-x_n|$ (la distanza tra le ultime 2 iterazioni).

Quando si ha a che fare con un metodo convergente, lo step va a zero per $n\to\infty$ perche' $x_n,x_{n+1}\to\xi$; ma una stima dell'errore basata sullo step non e' affidabile in generale (anche se  talvolta viene usata come stima empirica tipo "ultima spiaggia").

Invece col metodo di Newton lo step e' una stima molto buona per costruzione, almeno per $n$ abbastanza grande perche'
$$\text{step}(n)=|x_{n+1}-x_n|=\dfrac{f(x_n)}{f'(x_n)}$$
visto che $x_{n+1}=x_n-\dfrac{f(x_n)}{f'(x_n)}$, cioe' per Newton lo *step* e' per costruzione un **residuo pesato** ed e' gia' in pratica calcolato ad ogni iterazione visto che basta prendere il modulo di $\dfrac{f(x_n)}{f'(x_n)}$ che e' una quantita' chiave nel processo di calcolo.

In effetti, la stima a posteriori dello step e' generalmente affidabile, soprattutto quando e' accompagnata da altri controlli di convergenza.
Infatti, sia nel teorema di convergenza globale che in quello di convergenza locale si vede che *l'errore e' decrescente* ($e_{n+1}\lt e_n$): quindi per avere garanzie di essere in condizioni di convergenza si puo' controllare per prima cosa che lo step sia decrescente e poi usarlo come stima accurata dell'errore non appena $|f'(x_n)|$ comincia a stabilizzarsi, cioe' $\dfrac{f'(x_{n_1})}{f'(x_n)}\approx1$ (criterio empirico che abbiamo gia' discusso nella stima del residuo pesato per il metodo di bisezione).

Questo approccio puo' essere parcolarmente utile se si usa un metodo cosiddetto *ibrido*, come sono di solito i metodi adottati dai solutori automatici di di equazioni, che in input richiedono $f$, $x_0$ e una tolleranza, accoppiando un metodo "lento" ma affidabile (ad esempio la bisezione) per fare alcune iterazioni, partendo poi con un metodo veloce (come ad esempio Newton o le sue varianti) con un controllo di convergenza e iterando il procedimento finche' il metodo veloce non entra in condizioni di convergenza rapida perche' inizia lell'intorno "giusto" dell zero, arrestando le iterazioni con una stima a posteriori.

Sempre parlando di approccio empirico all'implementazione di un metodo iterativo convergente (non necessariamente un metodo per la soluzione di equazioni) siamo interesssati in pratica all'errore relativo, a meno che $\xi=0$ (o sia piccolissimo).

Per avere un test di arresto che funzioni sia per $\xi\neq0$ (dove conta l'errore relativo) sia per $\xi \approx 0$ (dove conta l'errore assoluto), avendo una stima affidabile dell'errore assoluto si usa spesso un test del tipo
$$\text{stima}(n)\le \epsilon_a+|x_n|\epsilon_r$$
dove $\epsilon_a$ e' una tolleranza assoluta ed $\epsilon_r$ una tolleranza relativa.
In questo modo, visto che $|x_n|\to|\xi|$, se $|\xi|$ e' ben staccato da zero ci si ferma su $\epsilon_r$ se $\xi=0$ (oppure $|\xi|$ e' "piccolissimo") ci si ferma su $\epsilon_a$.

Vediamo come si comporta la stima dello step per il metodo di Newton nel calcolo di $\sqrt2$:

![[img025.png]]

Si vede che la curva dello step e' *"parallela"* alla curva dell'errore (quindi e' un'ottima stima) ma risulta shiftata in avanti di $1$: questo e' naturale, visto che per avere il residuo pesato $\dfrac{|f(x_n)|}{|f'(x_n)|}$ al passo $n$ bisogna essere arrivati al passo $n+1$, $x_{n+1}-x_n=-\dfrac{f(x_n)}{f'(x_n)}$.
Siccome l'errore decresce velocemente, la stima dello step diventa ancora piu' affidabile
$$r_{n+1}=\dfrac{e_{n+1}}{|\xi|}\ll r_n = \dfrac{e_n}{|\xi|}\approx \dfrac{|x_{n+1}-x_n|}{|\xi|}$$

Facciamo ora due esempi dei applicabilita' del metodo di Newton.

### Esempio 1 (metodo di Erone per le radici quadrate)

Abbiamo visto come usare il metodo di Newton per calcolare $\sqrt2$ come soluzione dell'equazione algebrica (zero di un polinomio) $f(x) = x^2 - 2 =0$.
L'approccio e' generalizzabile al calcolo di $\sqrt{a}, a\gt0$, risolvendo l'equazione $f(x)=x^2-a=0$.

Osserviamo che $f(0)=-a\lt0$ e $f(b)\lt0$ per $b^2\gt a$.
Visto che $f'(x)=2x$ e inoltre $f''(x)=2\gt0$ siamo nelle ipotesi del teorema di convergenza globale scegliendo $x_0: {x_0}^2\gt a$.
Qual'e' la forma delle iterazioni di Newton in questo caso?
$$\begin{array}{c}
x_{n+1} = x_n - \dfrac{f(x_n)}{f'(x_n)}= x_n - \dfrac{{x_n}^2}{2x_n} \\
= \dfrac{2{x_n}^2 - {x_n}^2 + a}{2x_n} = \dfrac{{x_n}^2+a}{2x_n}\\
=\dfrac{x_n}{2} + \dfrac{a}{2x_n},\; n\ge0
\end{array}$$
Quindi se $a$ e' razionale (in particolare intero) e $x_0$ e' razionale il metodo fornisce per costruzione  una successione di razionali (frazioni) che converge a $\sqrt{a}$ (quadraticamente perche' $\sqrt{a}$ e' zero semplice).

Questa iterazione era gia' nota in eta' ellenistica ed e' attribuita al matematico greco Erone di Alessandria (che vi era arrivato non col calcolo differenziale, ignoto all'epoca, ma con metodi geometrici).

### Esempio 2 (applicazione del metodo di Newton a un'equazione trascendente)

Il metodo di Newton e' applicabile (nelle giuste ipotesi) a qualsiasi equazione del tipo $f(x)=0$, quindi anche ad equazioni in cui $f$ non e' un polinomio o una funzione razionale (rapporto di polinomi), dette equazioni trascendenti (quelle polinomiali/razionali sono dette equazioni algebriche).

E' il caso di osservare che in realta' anche le equazioni algebriche richiedono metodi numerici: infatti e' noto (teoria di Galois) che gli zeri di polinomi di grado $\ge5$ non sono calcolabili tramite radicali e d'altra parte le stesse equazioni di secondo grado richiedono il calcolo di $\sqrt\Delta$ che va fatto con un metodo approssimato (abbiamo visto sopra come usare Newton che e' molto veloce per le radici quadrate).

Consideriamo l'equazione trascendente $f(x)=x-e^{-\alpha x}=0, \alpha\gt0$ che si puo' interpretare come intersezione di grafici (in effetti ha banalmente la forma di punto fisso $x=e^{-\alpha x}$ e la useremo anche nel prossimo capitolo).

Osserviamo che $f\in C^\infty(\mathbb{R}), f(0)=-1\lt0, f(1)=1-e^{-\alpha}\gt0$ quindi $\exists\;\xi\in(0,1):f(\xi)=0$.
Lo zero e' sicuramente unico (in $\mathbb{R}$), $f'(x)=1+\alpha e^{-\alpha x} \gt 0$ ($f$ strettamente crescente), inoltre $f''(x)=-\alpha^2e^{-\alpha x}\lt0$ ($f$ strettamente concava) e infine $f'(x)\ge f'(1)=1+\alpha e^{-\alpha}$ e $|f''(x)|=\alpha^2e^{-\alpha x}\le \alpha^2, x\in[0,1]$.

![[img026.png]]

Vale quindi
$$c=\dfrac{1}{2}\dfrac{M_2}{m_1}\le \dfrac{1}{2}\dfrac{\alpha^2}{1+\alpha e^{-\alpha}}$$
per $\alpha\le1$ abbiamo che $c|\xi|\lt\dfrac{1}{2}$ da cui otteniamo $r_{n+1}=\dfrac{e_{n+1}}{|\xi|}\lt\dfrac{1}{2}{r_n}^2$ e anche in questo caso ci sara' un raddoppio (almeno) delle cifre corrette ad ogni iterazione.

Cocnludiamo la lezione mostrando in breve altri 2 metodi classici per la soluzione numerica di quazioni non lineari, anch'essi basati su una forma di *linearizzazione iterativa*, il metodo delle **corde** e il metodo delle **secanti**.

![[img027.png]]

Entrambi i metodi corrispondono a sostituire l'equazione $f(x)=0$ con un'equazione lineare del tipo
$$f(x_n)+q_n(x-x_n)=0$$
dove nel metodo delle corde $q_n$ e' il coefficiente angolare della corda (segmento) per $(x_n, f(x_n))$ e $(b,f(b))$
$$q_n=\dfrac{f(b)-f(x_n)}{b-x_n}$$
mentre nel metodo delle secanti $q_n$ e' il coefficiente angolare della retta per $(x_{n+1}, f(x_{n-1}))$ e $(x_n, f(x_n))$
$$q_n = \dfrac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}$$
(osserviamo che nel metodo di Newton $q_n=f'(x_n)$).

Entrambi i metodi hanno ipotesi di convergenza globale e locale, ad esempio il metodo delle corde converge se $f''(x)$ ha segno costante e $x_0$ e' tale che $f(x_0)f''(x_0)\gt0$; inoltre, il metodo delle corde ha convergenza lineare ($p=1$) mentre il metodo delle secanti ha convergenza superlineare con $p\in(1,2)$.
In effetti ci aspettiamo che il metodo delle secanti sia piu' veloce, perche' entrambi rispetto a Newton hanno un rapporto incrementale al posto della derivata
$$x_{n+1}=x_n-\frac{f(x_n)}{q_n}$$
con $n\ge0$ (corde) e $n\ge1$ (secanti).
Ma mentre nelle corde un estremo e' fisso, nelle secanti per $f\in C^1$
$$q_n=\frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}=f'(u_n)$$
(per il teorema del valor medio) con $u_n\in \text{int}(x_n, x_{n-1})$, quindi se c'e' convergenza per il teorema dei 2 carabinieri $u_n\to\xi, n\to\infty$ e $q_n\to f'(\xi), n\to\infty$, cioe' la secante tende ad essere sempre piu' "simile" ad una tangente al crescere di $n$.

Il metodo delle secanti risulta quindi una valida alternativa al metodo di Newton quando $f'$ non e' nota o difficile da calcolare.