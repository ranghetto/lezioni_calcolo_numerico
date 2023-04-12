## 10. Metodo di Newton: convergenza locale, ordine di convergenza, test di arresto, esempi, altri metodi di linearizzazione

Partiamo dalla relazione chiave ottenuta nel capitolo precedente per l'errore del metodo di Newton:
$$
\begin{array}{c}
e_{n+1}=c_ne_n^2 \leq ce_n^2 \\
c_n = \dfrac{1}{2}\dfrac{|f''(z_n)|}{|f'(x_n)|},\quad c=\dfrac{1}{2}\dfrac{M_2}{m_1} \\
M_2 = \underset{x\in[c,d]}{\max} |f''(x)|, \; m_1 = \underset{x\in[c,d]}{\min} |f'(x)|\gt0
\end{array}
$$