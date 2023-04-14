## 11. Punto Fisso: iterazioni punto fisso, teorema delle contrazioni, convergenza locale, ordine di convergenza, Newton come iterazione di punto fisso

In questa lezione studieremo equazioni della forma
$$x=\phi(x), x\in I\subseteq \mathbb{R}$$
dove $\phi\in C(I)$ e $I$ e' un intervallo chiuso (non necessariamente limitato) di $\mathbb{R}$ e la loro soluzione numerica tramite semplici iterazioni del tipo
$$x_{n+1}=\phi(x_n), n\ge0, x_0\in I$$
comunemente dette *iterazioni di punto fisso*.
In particolare, vedremo ipotesi che garantiscono la convergenza $x_n\to\xi, n\to\infty$ dove $\xi=\phi(\xi)$ e' detto punto fisso di $\phi\in I$ e studieremo l'ordine di convergenza dell'iterazione, scoprendo che il metodo di Newton puo' essere interpretato come iterazione di punto fisso.

Enunciamo qui sotto un famoso teorema, detto *teorema delle contrazioni*.

> **Teorema** (esistenza e unicita' del punto fisso e convergenza delle iterazioni di punto fisso per una _contrazione_)
> Sia $\phi:I\subseteq \mathbb{R}\to \mathbb{R}$ una funzione derivabile nell'intervallo chiuso $I\subseteq\mathbb{R}$, tale che:
> 1. $\phi(I)\subseteq I$ cioe' l'immagine di $I$ tramite $\phi$, $\phi(I)=\{y:y=\phi(x),x\in I\}$, e' contenuta in I;
> 2. $\exists\;\Theta\in(0,1):|\phi'(x)|\lt\Theta\;\forall\;x\in I \implies \exists!\xi\in I: \xi = \phi(\xi)$ (punto fisso) e $\forall\;x_0\in I, \xi= lim_{n\to\infty}x_n$ dove $x_{n+1}=\phi(x_n),n\ge0$ .

Prima di dimostrare questo teorema (che puo' essere esteso ad ambiti molto piu' astratti, qui ci limitiamo a funzioni reali di variabile reale, come provato dal matematico polacco S. Banach nel 1919 facneodlo diventare uno dei risultati chiave dell'analisi matematica contemporanea), facciamo alcune osservazioni:

1. l'intervallo $I$ e' assunto chiuso, ma puo' non essere limitato, cioe' $I=[a,b]$ con $-\infty\lt a\lt b \lt +\infty$ ma anche $I=[a,+\infty)$ oppure $I=(-\infty,b]$ o addirittura $I= \mathbb{R}$;
2. $\phi$ e' una contrazione (di $I$ in se' stesso), cioe' contrae le distanze di un fattore $\Theta\lt1$. Infatti per il teorema del valor medio $\forall\;x,y\in I$ vale $\phi(x)-\phi(y)=\phi'(z)(x-y), z\in \text{int}(x,y)\subset I$ da cui $|\phi(x)-\phi(y)|=|\phi'(z)||x-y| \le \Theta|x-y|\lt|x-y|$;
3. chiaramente la disuguaglianza appena provata, assunta direttamente come ipotesi, implica che $\phi$ e' continua in $I$, infatti $\forall\;x,\bar{x}\in I$ si ha che $0\le |\phi(x)-\phi(\bar{x})|\le\Theta|x-\bar{x}|$ e quindi per il teorema dei due carabinieri $|\phi(x)-\phi(\bar{x})|\to0,x\to\bar{x}$ che e' equivalente a dire che $\phi(x)\to\phi(\bar{x}), x\to\bar{x}$.

> **Dimostrazione** (del teorema)
> Cominciamo dimostrando l'esistenza di un punto fisso, limitandoci al caso di $[a,b]$ limitato:
> in questo caso basta l'ipotesi *1.* e la continuita' di $\phi$ (non server che $\phi$ sia una contrazione).
> Siccome $\phi$ e' continua, tale e' $f(x)=x-\phi(x)$ se $a=\phi(a)$ oppure $b=\phi(b)$ allora $a$ oppure $b$ sono punto fisso.
> Se invece $a\neq\phi(a)$ e $b\neq\phi(b)$ siccome $a\le\phi(x)\le b \;\forall\;x\in[a,b]$ si ha $a-\phi(a)\lt0$ e $b-\phi(b)\gt0$ cioe' $f$ e' continua e cambia segno agli estremi $\implies \exists\;\xi \in(a,b):f(\xi)=0$ cioe' $\exists\;\xi\in(a,b): \xi=\phi(\xi)$.
> La continuita' non basta pero' a garantire l'unicita' del punto fisso, come si vede da questo disegno
> ![[img028.png]]
>
>Ma se $\phi$ e' una contrazione, l'unicita' e' assicurata.
>Infatti se $\exists\xi_1,\xi_2\in I$ con $\xi_1\neq\xi_2$ tali che $\xi_1=\phi(\xi_1)$ e $\xi_2=\phi(\xi_2)$ allora $|\xi_1 - \xi_2| = |\phi(\xi_1)-\phi(\xi_2)|\le\Theta|\xi_1-\xi_2|$ cioe' $\Theta\ge1$ contro l'ipotesi che $\Theta\lt1$.
>Resta da dimostrare che $\forall\;x_0\in I$, definendo $x_{x+1}=\phi(x_n), n\ge0$ si ha $lim_{n\to\infty}x_n=\xi$ ora $e_{n+1}=|x_{n+1}-\xi|=|\phi(x_n)-\phi(\xi)|\le\Theta|x_n-\xi|=\Theta e_n$ da cui $e_1\le\Theta e_0, e_2\le\Theta e_1 \le \Theta^2 e_0, \dots, e_n\le \Theta^n e_0 \to 0, n\to\infty$ perche' $\Theta\in(0,1)$.

E' il caso di fare subito alcune osservazioni importanti:

1. nel teorema delle contrazioni, la dimostrazione generale e' basata sul fatto che la successione $\{x_n\}$ e' di Cauchy e quindi convergente a uno $\xi\in I$ (perche' $I$ e' chiuso) che 'e automaticamente punto fisso perche' per continuita' di $\phi$, $$\xi=\lim_{n\to\infty}x_{n+1}=\lim_{n\to\infty}\phi(x_n)=\phi(lim_{n\to\infty}x_n)=\phi(\xi)$$ Ma anche con la dimostrazione scritta sopra, si ottiene la **stima a priori** dell'errore $e_n\le\Theta^ne_0$. Se $\Theta$ (ed $e_0$) sono noti questa permette a priori di stabilire il numero di iterazioni sufficiente ad ottenere $\xi$ con una tolleranza $\epsilon\gt0$, risolvendo la disuguaglianza $$\Theta^ne_0\le\epsilon\iff e^{nlog\Theta}\le e^{log(\epsilon/e_0)}\iff n\ge \dfrac{\log(\epsilon/e_0)}{\log(\Theta)}=-\dfrac{\log(\epsilon/e_0)}{|\log\Theta|}=\dfrac{\log(e_0/\epsilon)}{|\log\Theta|}$$ dato che $\Theta\in(0,1)$ e $\log\Theta\lt0$. Notiamo che la convergenza e la stima ottenute valgono $\forall x_0\in I$, cioe' le iterazioni di punto fisso che costruiscono (infinite) successioni diverse l'una dalll'altra al variare di $x_0$, in ogni caso forniscono successioni che convergono tutte allo *stesso limite* che e' *l'unico punto fisso* di $\phi$ in $I$, con una convergenza che e' almento lineare perche' $e_{n+1}\le\Theta e_n$;
2. si puo' facilmente ottenre una **stima a posteriori** dell'errore che spesso e' piu' precisa della stima a priori: basta infatti scrivere $$x_{n+1}-\xi=x_{n+1}-x_n+x_n-\xi$$ ma $$x_{n+1}-\xi=\phi(x_n)-\phi(\xi) = \phi'(z_n)(x_n -\xi),\quad z_n\in \text{int}(x_n,\xi)$$ da cui $$\phi'(z_n)(x_n -\xi) = x_{n+1}-x_n+x_n-\xi$$ ovvero $$\begin{array}{c}(1-\phi'(z_n))(x_n-\xi)=x_n-x_{n+1}\\ \implies \dfrac{|x_{n+1}-x_n|}{1-\phi'(z_n)}=e_n\le\underbrace{\dfrac{|x_{n+1}-x_n|}{1-\Theta}}_{\text{prima stima a posteriori}}\end{array}$$ perche' $|\phi'(z_n)|\le\Theta\lt1$ $$\implies |1-\phi'(z_n)|\ge|1-|\phi'(z_n)||\ge1-\Theta\gt0\quad \text{e} \quad \dfrac{1}{|1-\phi'(z_n)|}\le\dfrac{1}{1-\Theta}$$
In pratica abbiamo fatto vedere che l'errore e' stimato dallo *step* $=|x_{n+1}-x_n|$ a meno del fattore $\dfrac{1}{1-\Theta} =$ *peso*.
Se $\Theta$ e' piccolo lo step diventa da solo una buona stima dell'errore perche' il perso e' $\approx 1$; invece se $\Theta$ e' vicino ad $1$ lo step va corretto per evitare una possibile sottostima, ma se $\phi\in C^1(I)$, siccome $z_n\to\xi$ allora $\phi'(z_n)\to\phi'(\xi), n\to\infty$, quindi una stima a posteriori migliore e' tendenzialmente la stima empirica (almeno per $n$ abbastanza grande) 
$$e_n=\dfrac{|x_{n+1}-x_n|}{1-\phi'(z_n)}\approx\underbrace{\dfrac{x_{n+1}-x_n}{1-\phi'(x_n)}}_{\text{seconda stima a posteriori}}$$

### Esempio
Consideriamo l'equazione di punto fisso
$$x=\phi(x)=e^{-\alpha x}, \alpha\gt0$$
se $\alpha \lt 1, |\phi'(x)|=|-\alpha e^{-\alpha x}|\le\alpha\lt1$ e quindi $\phi$ contrae le distanze, quindi l'ipotesi *2.* del teorema delle contrazioni e' soddisfatta con $\Theta=\alpha$; d'altra parte, $0\lt\phi(x)\le1, \forall x \in [0,+\infty)$, quindi anche l'ipotesi $1.$ e' soddisfatta con $I=[0,+\infty)$.
In realta' siccome $\phi'(x)=-\alpha e^{-\alpha x}\lt0$, $\phi$ e' strettamente descrescente (e positiva), quinid visto che $\phi(0)=1$ e $0\lt\phi(1)=1-e^{-\alpha}\lt1$ si ha $0\lt\phi(x)\le1, \forall x\in[0,1]$ cioe' $\phi$ e' anche una contrazione di $[a,b]=[0,1]$ in se stesso.

Allora $\exists!\;\xi$ punto fisso di $\phi\in[0,1]$: prendiamo $x_0 = \dfrac{1}{2} \implies e_0=|x_0-\xi|\le\dfrac{1}{2}$ e la successione $x_{n+1}=\phi(x_n),\;n\ge0$ converge a $\xi$ con la stima a priori dell'errore $e_n\le\dfrac{\alpha^n}{2}$.
Nel grafico qui sotto mostriamo l'errore effettivo, la stima a priori e la stima a posteriori dello step (non pesato) e dello step pesato da $\dfrac{1}{1-\alpha}$ e da $\dfrac{1}{1-\phi'(x_n)}$ con
$$\begin{array}{c}
\alpha=0.2\implies \text{fl}(\xi)=0.8445798\dots\\
\alpha=0.9\implies \text{fl}(\xi)=0.5887032\dots
\end{array}$$
(questi valori sono stati calcolati col metodo di Newton)

![[img029.png]]

Si vede chiaramente che la convergenza e' lineare, che la stima a priori e' una sovrastima, molto distante dall'errore per $\alpha=0.9$: infatti $\dfrac{e_{n+1}}{e_n}=|\phi'(z_n)|\to|\phi'(\xi)|,n\to\infty$ cioe' $|\phi'(\xi)|=\alpha e^{-\alpha \xi}=L$ e' la costante asintotica, il parametro che effettivamente regola la velocita' di convergenza ($\alpha$ e' solo una stima) perche' per $n$ abbastanza grande $e_{n+1}\approx Le_n$.
Qui per $\alpha=0.9$ si ha $L\approx0.53$ che e' ben minore di $\alpha$, mentre per $\alpha=0.2$ si ha $L\approx0.17$ che e' poco minore di $\alpha$.

In effetti la seconda stima a posteriori, $e_n\approx\dfrac{|x_{n+1}-x_n|}{1-\phi'(x_n)}$, e' una stima aderente dell'errore ma e' shiftata in avanti di $1$, fenomeno che abbiamo gia' visto nella stima con lo step nel motodo di Newton, perche' per stimare $e_n$ bisogna essere al passo $n+1$ (*step=$|x_{n+1}-x_n|$*).

Dopo aver discusso questo esempio semplice ma significativo, vediamo che anche per le iterazioni punto fisso vale un risultato di *convergenza locale* (mentre la formulazione generale del teorema delle contrazioni ha carattere *globale*, visto che $x_0\in I$ e' arbitrario e per la convergenza non e' importante che $x_0$ sia vicino a $\xi$).

> **Teorema** (convergenza locale delle iterazioni di punto fisso)
> Sia $\xi$ punto fisso di $\phi\in C^1(I_\delta(\xi))$ dove $I_\delta(\xi)=[\xi-\delta, \xi+\delta], \delta \gt0$ e sia $|\phi'(\xi)|\lt1$ allora $$\exists \delta' \le \delta:x_{n+1}=\phi(x_n), n\ge0$$ converge a $\xi, \forall x_0 \in I_{\delta'}(\xi)$.

E' chiaro il carattere *locale* di questo risultato, che fornisce condizioni sufficienti per la convergenza delle iterazioni di punto fisso purche' $x_0$ *sia abbastanza vicino a* $\xi$.

Come per tutti i metodi iterativi, e' importante capire quale sia l'*ordine di convergenza* delle iterazioni di punto fisso.
Abbiamo gia' osservato che nel caso di una contrazione l'ordine e' almento $p=1$ perche' vale $e_{n+1}\le\Theta e_n$ con $\Theta\in(0,1)$.
D'altra parte, nell'esempio svolto prima abbiamo fatto vedere che l'ordine e' esattamente $p=1$ se $\phi'(\xi)\neq0$ con costante asintotica $L=|\phi'(\xi)|$. Diamo ora una caratterizzazione completa col seguente teorema.

>**Teorema** (ordine di convergenza delle iterazioni di punto fisso)
>Sia $\xi$ punto fisso di $\phi\in C^p(I),p\ge1$, con $I$ intervallo di $\mathbb{R}$ e supponiamo di essere in ipotesi che garantiscano la convergenza a $\xi$ di $x_{n+1}=\phi(x_n), n\ge0, x_0\in I$ (ad esempio le ipotesi del teorema delle contrazioni). Allora:
>1. $\{x_n\}$ ha ordine esattamente $p=1 \iff 0\lt|\phi'(\xi)|\lt1$;
>2. $\{x_n\}$ ha ordine esattamente $p\gt1 \iff \phi^{(j)}(\xi)=0, 1\le j\le p-1$ e $\phi^{(p)}(\xi)\neq0$.

>**Dimostrazione**
>$1.$ si dimostra subito visto che $e_{n+1}=|\phi'(z_n)|e_n, z_n\in\text{int}(\xi, x_n)$ per il teorema del valor medio, quindi $\lim_{n\to\infty}\dfrac{e_{n+1}}{e_n}=|\phi'(\lim z_)|=|\phi'(\xi)|$.
>Per il $2.$ utilizziamo la formula di Taylor di grado $p-1$ centrata in $\xi$ con resto $p$-esimo in forma di Lagrange $$\begin{array}{c}x_{n+1}=\phi(x_n)=\phi(\xi)+\phi'(\xi)(x_n-\xi)+\dfrac{\phi''(\xi)}{2}(x_n-\xi)^2+\cdots\\\cdots+\dfrac{\phi^{(p-1)}(\xi)}{(p-1)!}(x_n-\xi)^{p-1}+\dfrac{\phi^{(p)}(u_n)}{p!}(x_n-\xi)^p, u_n\int\text{int}(\xi,x_n)\end{array}$$ Dimostriamo prima "$\Longleftarrow$" (condizione sufficiente): se $\phi^{(j)}(\xi)=0, 1\le j\le p-1$ e $\phi^{(p)}(\xi)\neq0$, da Taylor $x_{n+1}-\xi=\dfrac{\phi^{(p)}(u_n)}{p!}(x_n-\xi)^p$ e passando ai moduli $$\dfrac{e_{n+1}}{e_n^p}=\dfrac{|\phi^{(p)}(u_n)|}{p!}\to\dfrac{\phi^{(p)}(\xi)}{p!}\neq0, n\to\infty$$ perche' $u_n\to\xi, n\to\infty$ e $\phi^{(p)}$ e' continua quindi $\lim\phi^{(p)}(u_n) = \phi^{(p)}(\lim u_n) = \phi^{(p)}(\xi)$ ovvero $\{x_n\}$ ha ordine esattamente $p$.
>Per dimostrare "$\Longrightarrow$" (condizione necessaria) supponiamo per assurdo che $\{x_n\}$ abbia ordine esattamente $p$ ma che $\exists j\lt p$ tale che $\phi^{(j)}\neq0$.
>Prendiamo $k=\min\{j\lt p:\phi^{(j)}(\xi)\neq0\}$ e scriviamo $\dfrac{e_{n+1}}{e_n^p}=\dfrac{e_{n+1}}{e_n^k}e_n^{k-p}$. Ora per ipotesi $\dfrac{e_{n+1}}{e_n^p}\to L\neq0$, d'altra parte, con lo stesso ragionamento usato per la dimostrazione della condizione sufficiente, tramite la formula di Taylor si avrebbe $$\dfrac{e_{n+1}}{e_n^k}\to\dfrac{|\phi^{(k)}(\xi)|}{k!}=L'\neq0$$ ma allora $\dfrac{e_{n+1}}{e_n^p}=\underbrace{\dfrac{e_{n+1}}{e_n^k}}_{L'}\underbrace{e_n^{k-p}}_{\infty}$ perche' $k-p\lt0$ e $e_n\to0$ cioe' alla fine $\dfrac{e_{n+1}}{e_n^p}\to\infty, n\to\infty$, contraddicendo l'ipotesi che abbia limite finito.

Ribadiamo che la condizione data e' *necessaria* e *sufficiente*, cioe' fornisce una caratterizzazione completa di quando le iterazioni di punto fisso hanno ordine $p\ge1$.
In particolare (e questo ci servira' tra poco) le iterazioni di punto fisso possono avere convergenza *quadratica*, *cubica*, etc...

Per concludere il capitolo mostriamo infatti che il metodo di Newton si puo' re-interpretare come iterazione di punto fisso e che in base a quanto visto sopra si vede subito che ammette convergenza locale e che la convergenza e' quadratica per zeri semplici.
Infatti l'iterazione di Newton $$x_{n+1}=x_n-\dfrac{f(x_n)}{f'(x_n)},n\ge0$$ e' di tipo punto fisso con $\phi(x)=x-\dfrac{f(x)}{f'(x)}$ e se $f\in C^k(I)$ con $f'(x)\neq0, \forall x\in I$ intervallo, allora $\phi\in C^{k-1}(I)$.
E' evidente che $\xi$ e' zero (semplice) di $f\iff\xi$ e' punto fisso di $\phi$.
Inoltre, il teorema di convergenza locale per Newton (zeri semplici) discende immediatamente dal teorema di convergenza locale per le iterazioni di punto fisso se $f\in C^2(I_\delta(\xi))$, infatti
$$\phi'(x)=\dfrac{d}{dx}(x-\dfrac{f}{f'})=1-\underbrace{\frac{(f')^2-ff''}{(f')^2}}_{\text{derivata del rapporto}}=\dfrac{ff''}{(f')^2}$$quindi $f(\xi)=0\implies\phi'(\xi)=0$ e $|\phi'(\xi)|=0\lt1$, allora $\exists\delta'\le\delta$ tale che l'iterazione di Newton converge come iterazione di punto fisso $\forall x_0\in I_{\delta'}(\xi)$.

D'altra parte, sempre interfpretando Newton come iterazione di punto fisso, e' immediato che la convergenza per zeri semplici e' almeno quadratica perche' $\phi'(\xi)=0$ ed e' esattamente quadratica se $f''(\xi)\neq0$, utilizzando la caratterizzazione vista sopra (che pero' richiede $\phi\in C^2$ e quindi $f\in C^3$).
Infatti $$\phi''=\dfrac{d}{dx}(\dfrac{ff''}{(f')^2})=\dfrac{(f'f''+ff''')(f')^2-2f'f''(ff'')}{(f')^4}$$ da cui $\phi''(\xi)=\dfrac{f''(\xi)}{f'(\xi)}\neq0$ perche' $f''(\xi)\neq0$ (e $f'(\xi)\neq0$).

Non e' difficile vedere che interpretando Newton come iterazione di punto fisso, nel caso di zero multiplo l'ordine di convergenza diventa $p=1$ con costante asintotica $|\phi'(\xi)|=1-\dfrac{1}{m}$ dove $m$ e' la molteplicita' di $\xi$ (il numero di derivate successive che si annullano in $\xi$) e che se $m$ e' nota l'iterazione
$$x_{n+1}=x_n-m\dfrac{f(x_n)}{f'(x_n)}=\phi_m(x_n)$$ torna di ordine almeno $p=2$ perche' $\phi'_m=\lim_{x\to\xi}\phi'_m(x)=0$.