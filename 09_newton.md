## 9. Metodo di Newton (tangenti) convergenza e velocita' di convergenza

Come abbiamo visto nello scorso capitolo, il metodo di bisezione e' un metodo semplice ed efficace per la soluzione numerica di equazioni non lineari, in grado di funzionare con richieste minimali, sia dal punto di vista analitico (con l'ipotesi del teorema degli zeri delle funzioni continue) sia computazionali (basta saper calcolare il segno di $f(x_n$) cioe' fare su $f(x_n)$ un errore relativo $\lt 100\%$).

Accanto a questi indubbi vantaggi, il metodo di bisezione ha pero' un handicap: e' abbastanza "lento".
Come abbiamo visto, la stima a priori dell'errore decade di un fattore $\dfrac{1}{2}$ ad ogni iterazione e l'errore lo fa "in media", cioe'
$$e_{n+1} \approx \dfrac{1}{2}e_n$$
non esattamente ad ogni iterazione ma in media su un certo numero di iterazioni.

Ovviamente la stessa proprieta' vale in media per l'errore relativo
$$r_{n+1} = \dfrac{e_{n+1}}{|\xi|}\approx\dfrac{1}{2}\dfrac{e_n}{|\xi|}\dfrac{1}{2}r_n$$
per $\xi\neq0$; sono quindi necessarie in media 3-4 iterazioni affinche' l'errore relativo scenda di $\dfrac{1}{10}$ (cioe' per guadagnare una cifra decimale corretta), visto che
$$\dfrac{1}{16} = \Big(\dfrac{1}{2}\Big)^4 \lt \dfrac{1}{10} \lt \Big(\dfrac{1}{2}\Big)^3 = \dfrac{1}{8}$$
(infatti abbiamo visto che per calcolare $\sqrt{2}$ alla precisione di macchina in Matlab servono una cinquantina di iterazioni e in effetti $\epsilon_M=2^{-53} \approx 10^{-16}$).

In questa lezione introdurremo un metodo molto piu' efficiente per il calcolo di zeri, ovvero il *metodo di Newton*, o *metodo delle tangenti*, che come si capisce dal nome risale agli albori del calcolo differenziale nel $XIIV$ secolo.

La maggiore efficienza avra' un prezzo in termini di richieste analitiche ($f$ almeno derivabile in $[a,b]$) e computazionali (saper calcolare sia $f$ sia $f'$ con buona accuratezza).

L'idea del metodo e' semplice: si tratta di **linearizzare** iterativamente l'equazione $f(x)=0$ sostituendo $f$ ad ogni iterazione con la retta tangente nel punto $(x_n, f(x_n))$ del grafico (purche' $f$ sia derivabile), come vediamo nel disegno qui sotto

![[img018.png]]

Per trovare l'espressione analitica delle iterazioni, calcoliamo lo zero della retta tangente nel punto $(x_n, f(x_0))$:
$$\begin{cases}
y=0 && \text{asse }x \\
y=f(x_0) + f'(x_0)(x-x_0) && \text{retta tangente}
\end{cases}$$
dove stiamo usando l'interpretazione geometrica della derivata nel punto $x_0$ come coefficiente angolare della retta tangente nel corrispondente punto del grafico di $f$.
Otteniamo l'equazione $o = f(x_0)+f'(x_0)(x-x_0)$ la cui soluzione e'
$$x_1 = x_0 - \dfrac{f(x_0)}{f'(x_0)},\quad f'(x_0)\neq0$$
Osserviamo subito che se $f'(x_0)=0$ la retta tangente sarebbe parallela all'asse $x$ e non avrebbe un punto di intersezione connesso; ma e' essenziale anche la posizione del valore iniziale $x_0$, se scelto "male" gia' la prima iterazione potrebbe far uscire dall'intervallo di definizione:

![[img019.png]]

In generale, per ottenere $x_{n+1}$ a partire da $x_n$ (se $f'(x_n)\neq0$) si cerca l'intersezione della tangente in $(x_n, f(x_n))$ con l'asse $x$ ($y=0$)
$$\begin{cases}
y=0 \\
y=f(x_n)+f'(x_n)(x-x_n)
\end{cases}$$
ottenendo la formula iterativa
$$x_{n+1} = x_n - \dfrac{f(x_n)}{f'(x_n)},\quad n=0,1,2,\dots$$
Ora, come per tutti i metodi che producono una successione, cerchiamo delle condizioni che garantiscano la *convergenza* (in questo caso il limite dovra' essere lo zero $\xi$ di $f$).

Ci sono vari set di condizioni *sufficienti* che garantiscono la convergenza del metodo di Newton; ne mostriamo uno che assume altre due ipotesi del teorema degli zeri:
- la derivabilita', essenziale perche' la curva grafico ammetta tangente in ogni punto;
- convessita' o concavita' stretta, tramite il segno di $f''$.

> **Teorema** (convergenza del metodo di Newton con $f''$ di segno costante)
> Sia $f:[a,b]\to\mathbb{R}$ derivabile due volte in $[a,b]$, $f(a)f(b)\lt0$, $f''(x)\gt0\;\forall\;x\in[a,b]$ oppure $f''(x)\lt0\;\forall\;x\in[a,b]$, $x_0\in[a,b]$ tale che $f(x_0)f''(x_0)\gt0$ allora il metodo di Newton e' ben definito (cioe' $f'(x_n)\neq0\;\forall\;n$) e converge all'unico zero $\xi$ di $f$ in $(a,b)$.

 Ci sono 4 casi possibili in base al segno di $f''$ ovvero
 ![[img020.png]]

in $1)$ e $2)$ $f$ e' strettamente convessa, in $3)$ e $4)$ concava, in $1)$ e $3)$ $x_0$ va scelto in $(\xi, b]$, mentre in $2)$ e $4)$ $x_0$ va scelto in $[a, \xi)$.
Vediamo dai disegni che non e' escluso che $f'$ possa cambiare segno, l'ipotesi chiave e' che non cambi segno $f''$; ovviamente sono compresi i casi in cui $f''$ non cambia segno e anche $f'$ non cambia segno, cioe' ad esempio quando $f$ sstrettamente convessa e strettamente crescente.
Per semplicita' trattiamo solo il caso $1)$ perche' negli altri casi la dimostrazione e' analoga

>**Dimostrazione** (caso $1)$)
>Siamo quindi in questa situazione
>![[img021.png]]
>
>con $f(a)\lt0, f(b)\gt0, f''(x)\gt0\;\forall\;x\in[a,b], x_0\in(\xi,b]$.
>La dimostrazione si fa per induzione, mostrando che se $x_n\in(\xi,b]$ anche $x_{n+1} \in(\xi,b]$ e inoltre $x_{n+1}\lt x_n$. Infatti $x_{n+1}=x_n-\dfrac{f(x_n)}{f'(x_n)}$.
>Ma se $x_n\in(\xi,b]$ allora $f'(x_n)\gt0$ e $f(x_n)\gt0$, quindi $x_{n+1}$ si ottiene da $x_n$ sottraendo una quantita' $\gt0$, e quindi $x_{n+1}\lt x_n$.
>D'altra parte $f$ e' strattamente convessa, il che e' equivalente a dire che la tangente sta "sotto al grafico" $\forall\;x\in[a,b]$.
>Ma allora la tangente in un punto $\in(\xi,b]$ interseca l'asse $x$ a destra di $\xi$, cioe' se $x_n\in(\xi,b]$ anche $x_{n+1}\in(\xi,b]$.
>
>In definitiva, abbiamo provato che la successione ${x_n}$ e' decrescente e che $s_n\gt\xi\;\forall\;n$.
>Dall'analisi matematica e' noto che una succesisone monotona e limitata ha limite e che il limite e' $\sup\{x_n\}$ se e' crescente e $\inf\{x_n\}$ se e' decrescente (che corrisponde al nostro caso).
>Quindi $\exists\;\lim_{n\to\infty}x_n = \inf\{x_n\} = \eta$ con $\eta\ge\xi$.
>Ricordiamo infatti che le disuguaglianze conservano il verso passando al limite oppure facendo $\sup$ e $\inf$.
>Per concludere la dimostrazione, basta passare al limite nella formula che definisce il metodo (per brevita' scriveremo $\lim$ per $\lim_{n\to\infty}$)
>
>$$\begin{array}{c}
\eta=\lim x_{n+1} = \lim x_n-\dfrac{f(x_n)}{f'(x_n)} = \\
= \lim x_n - \lim \dfrac{f(x_n)}{f'(x_n)}=\\
= \lim x_n - \dfrac{\lim f(x_n)}{\lim f'(x_n)}=\\
= \lim x_n - \dfrac{f(\lim x_n)}{f'(\lim x_n)}=\\
= \eta - \dfrac{f(\eta)}{f'(\eta)}
\end{array}$$
> dove abbiamo usato le proprieta' dei limiti e la continuita' di $f$ ed $f'$ (portando il limite "dentro le funzioni").
> Quindi $\eta = \eta - \dfrac{f(\eta)}{f'(\eta)}=0 \implies f(\eta)=0$
> Ma allora $\eta = \xi$, perche' nelle ipotesi fatti (teorema degli zeri e $f''$ di segno costante) lo zero e' unico, quindi il metodo di Newton e' ben definito e $\{x_n\}$ converge a $\xi$.

E' il caso di ribadire che quello che abbiamo utilizzato e' uno dei vari set di condizioni sufficienti (ce ne sono altri), con ipotesi tipicamente di tipo "geometrico" (qui segno costante di $f''$ quindi $f$ e' strettametne convessa o concava) che garantiscono una convergenza che possiamo chiamare "globale" (cioe' non e' importante quanto il punto iniziale $x_0$ sia vicino allo zero, purche' sia nella zona giusta dell'intervallo).

Vedremo che il metodo di Newton puo' convergere anche con condizioni meno forti, purche' $x_0$ sia scelto in un intorno opportuno di $\xi$ (convergenza che chiameremo "locale").
Quest'ultimo e' uno dei punti di forza del metodo, perche' ad esempio la richiesta fatta sopra che $f''$ abbia segno costante e' piuttosto restrittiva e limiterebbe fortemente la classe di equazioni risolvibili.

L'altro essenziale punto di forza del metodo di Newton e' la *velocita' di convergenza*.

Per cominciare ad apprezzare questo apetto (che poi studieremo in dettaglio), facciamo un esempio in cui confrontiamo il metodo di Newton col metodo di bisezione


###### Esempio (calcolo di $\sqrt{2}$ con bisezione e con Newton)

Abbiamo gia' studiato il calcolo di $\sqrt{2}$ alla precisione di macchina col metodo di bisezione, applicabile in $[a,b]=[1,2]$ dove sono soddisfatte le ipotesi del teorema degli zeri per $f(x)=x^2-2$.
Inoltre $f'(x)=2x$ e $f''(x)=2\gt0$ (si tratta di un ramo di parabila con concavita' verso l'alto) quindi il metodo di Newton e' applicabile, partendo ad esempio da $x_0=2$ (essendo soddisfatte tutte le ipotesi del teorema dimostato prima, in particolare $f(x_0)f'(x_0)\gt0$).

Mostriamo ora la sequenza di iterazioni della bisezione, sapendo che in doppia precisione $fl(\sqrt2)=1.414213562373095$ in cui riquadriamo le cifre corrette
$$\begin{array}{c}
x_0 = \fbox1.5, x_1=\fbox1.25, x_2=\fbox1.375, x_3=\fbox{1.4}375, x_4=\fbox{1.4}0625, x_5=\fbox{1.4}21875,\\ x_6=\fbox{1.41}40625, \dots, x_{50}=fl(\sqrt{2})
\end{array}$$
Come ci aspettiamo ci vogliono 3-4 iterazioni per  guadagnare una cifra decimale corretta e una cinquantina di iterazioni per averne 16 corrette.
Invece con Newton
$$\begin{array}{c}
x_0 = 2, x_1=\fbox1.5, x_2=\fbox{1.41}6\dots67, x_3=\fbox{1.41421}5686174510,\\
x_4=\fbox{1.41421356237}4690, x_5 = fl(\sqrt{2})
\end{array}$$
cioe' con 5 iterazioni si ottiene $\sqrt{2}$ alla precisione di macchina.

Si puo' notare come il numero di cifre decimali corrette stia sostanzialmente *raddoppiando* ad ogni iterazione. In effetti $x_2$ ne ha 3, $x_3$ ne ha 6, $x_4$ ne ha 12 e $x_5$ ne ha addirittura 24 (se usassimo una precisione estesa, cioe' oltre le 16 cifre decimali).

Questo mostra che il metodo di Newton puo' essere estremamente piu' veloce del moetodo di bisezione, che pure ha una convergenza di tipo esponenziale, con errore proporzionale a $\Big(\dfrac{1}{2}\Big)^n$, cioe' Newton puo' convergere *piu' che esponenzialmente*.
Dal punto di vista dei grafici di errore (in scala logaritmica) la situazione nel alcolo di $\sqrt{2}$ e' la seguente:

![[img024.png]]

Perche' il metodo di Newton e' cosi' veloce? Per capirlo, dobbiamo analizzare l'errore, in particolare quale relazione leghi $e_{n+1}$ con $e_n$ (dove come al solito $e_n=|x_n-\xi|$).

Enunciamo il risultato sul comportamento dell'errore come teorema, che poi dimostreremo.

> **Teorema** (sulla velocita' di convergenza del metodo di Newton)
> Sia $f\in C^2[a,b]$ e si assuma di essere in ipotesi che garantiscono la convergenza del emtodo di Newton a $\xi \in [a,b]:f(\xi)=0$.
> Sia inoltre $\{x_n\}\subset[c,d]\subseteq[a,b]$ con $f'(x)\neq0\;\forall\;x\in[c,d]$
> Queste ipotesi implicano 
> $$\implies e_{n+1}\leq c\cdot e^2_n, n\geq0, c=\dfrac{1}{2}\dfrac{M_2}{m_1}$$
> con $M_2 = \underset{x\in[c,d]}{\max}|f''(x)|, m_1=\underset{x\in[c,d]}{f'(x)}\gt0$

Prima di dimostrare questa disuguaglianza, osserviamo che:
1. l'ipotesi $\{x_n\}\subset[c,d]$ con $f'(x)\neq0$ in $[c,d]$, come abbiamo gia' visto nell'analisi del residuo pesato col metodo di bisezione, ci assicura che lo zero $\xi$ e' semplice, cioe' $f'(\xi)\neq0$;
2. tale ipotesi e' soddisfatta ad esempio nelle condizioni del teorema di convergenza dimostrato prima con $[c,d]=[\xi,b]$ nei casi $1)$ e $3)$, $[c,d]=[a,\xi]$ nei casi $2)$ e $4)$.

> **Dimostrazione**
> Applicando la formula di Taylor centrata in $x_n$ e calcolata in $\xi$, con resto del secondo ordine in forma di Lagrange otteniamo
> $$\underset{=0}{f(\xi)}=f(x_n)+f'(x_n)(\xi-x_n)+\dfrac{f''(z_n)}{2}(\xi - x_n)^2$$
> dove $z_n \in int(x_n, \xi)\subset[c,d]$, da cui
> $$-\dfrac{f(x_n)}{f'(x_n)}=\xi-x_n+\dfrac{f''(z_n)}{2f'(x_n)}(\xi-x_n)^2$$
> Ma dalla definizione del metodo $-\dfrac{f(x_n)}{f'(x_n)}=x_{n+1}-x_n$ che inserito nelal formula di Tayloer porta a 
> $$x_{n+1}-\cancel{x_n}=\xi-\cancel{x_n}+\dfrac{f''(z_n)}{2f'(x_n)}(\xi-x_n)^2$$
> ovvero mettendo i moduli $e_{n+1} = |x_{n+1}-\xi|=c_n\cdot e_n^2$ con $c_n= \dfrac{1}{2}\dfrac{|f''(z_n)|}{f'(x_n)}$.
> La successione $\{c_n\}$ e' limitata, infatti $|f''(z_n)| \leq \underset{x\in[c,d]}{\max}|f''(x)| = M_2$, applicando il teorema di Weierstrass sull'$\exists$ di $\max$ e $\min$ assoluti a $|f''(x)|\in C[c,d]$; d'altra parte, applicando lo steso teorema a $|f'(x)|\in C[c,d]$ abbiamo che $m_1 = \underset{x\in[c,d]}{\min}|f'(x)|\gt0$ (perche' $\exists\;\bar{x}: m_1 = |f'(\bar{x})|$ e $f'(\bar{x})\neq0$).
> Otteniamo quindi $|f'(x_n)|\geq m_1 \gt0$ e infine $c_n\leq\dfrac{1}{2}\dfrac{M_2}{m_1} = c$

La relazione di tipo *quadratico* $e_{n+1}\leq c\cdot e_n^2$ e' la chiave per spiegare la velocita' di convergenza del metodo di Newton, perche' ci dice in sostanza che l'errore al passo $n+1$ e' maggiorato da una quantita' *proporzionale* (con costante di porporzionalita' non dipendente da $n$) al *quadrato delll'errore al passo* $n$.
Si noti la notevole differenza col metodo di bisezione, dove la relazione tra $e_{n+1}$ ed $e_n$ e' lineare ($e_{n+1}\approx\dfrac{1}{2}e_n$, in media).

Per apprezzare l'effetto della relazione quadratica, prima di tutto osserviamo che
$$c\cdot e_{n+1}\leq c\cdot c \cdot e_n^2=(c\cdot e_n)^2$$
Ora, fissiamo $\Theta\in(0,1)$: siccome abbiamo assunto che il metodo sia convergente, $ce_n\to0$, $n\to\infty$ e quindi $\exists\;\bar{n}: ce_n\leq\Theta\;\forall\;n\geq\bar{n}$ (con $\bar{n}$ dipendente da $\Theta$).
Applicando la disuguaglianza $ce_{n+1}\leq(ce_n)^2$ per $n\geq\bar{n}$:
$$\begin{array}{c}
ce_{\bar{n}+1} \leq (ce_\bar{n})^2 \leq \Theta^2\\
ce_{\bar{n}+2} \leq (ce_\bar{n+1})^2 \leq (\Theta^2)^2 = \Theta^4 \\
ce_{\bar{n}+3} \leq (ce_\bar{n+2})^2 \leq (\Theta^4)^2 = \Theta^8 \\
\vdots \\
ce_{\bar{n}+k} \leq (ce_\bar{n+k-1})^2 \leq (\Theta^{2^{k-1}})^2 = \Theta^{2^k} \\
\end{array}$$
Adesso, solo per fissare le idee nel confronto col metodo di bisezione, prendiamo $\Theta =\dfrac{1}{2}$.
Otteniamo, dopo $k$ iterazioni di Newton a partire da $\bar{n}$
$$e_{\bar{n}+k}^{\text{Newton}} \leq \dfrac{1}{c}\Big(\dfrac{1}{2}\Big)^{2^k}$$
mentre con $k$ iterazioni del metodo di bisezione
$$e_k^{bisezione}\leq\approx\Big(\dfrac{1}{2}\Big)^ke_0^{\text{bisezione}}$$
Il confronto sulle $k$ iterazioni va fatto guardando gli esponenti di $\dfrac{1}{2}$: nella bisezione l'esponente e' $k$ (cioe' cresce linearmente in $k$), con Newton invece l'esponente e' $2^k$ (cioe' cresce *esponenzialmente* in $k$).
Ad esempio per $k=6$ nella stima per la bisezione compare $\Big(\dfrac{1}{2}\Big)^6 = \dfrac{1}{64}\approx1.6\cdot 10^{-2}$ mentre nella stima per Newton $\Big(\dfrac{1}{2}\Big)^{2^6}=\Big(\dfrac{1}{2}\Big)^{64}\approx5\cdot10^{-20}$
Ribadiamo che $\Theta=\dfrac{1}{2}$ e' stato scelto arbitrariamente solo per fare un confronto diretto col metodo di bisezione, dove il fattore di riduzione $\dfrac{1}{2}$ e' intrinseco nella costruzione; nel caso del metodo di Newton (per zeri semplici) possiamo dire sostanzialmente che non appena $ce_n\lt1$ si innesca una riduzione rapidissima dell'errore, di tipo *"quadrati successivi"* che viene detto **convergenza quadratica**.