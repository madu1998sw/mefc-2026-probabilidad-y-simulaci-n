# Memoria de Resolución
## Máster Executive en Finanzas Cuantitativas 2026 — AFI Global Education
### Fundamentos Matemáticos: Probabilidad y Simulación

---

> Este documento es la memoria explicativa del trabajo de grupo. Recoge el desarrollo teórico completo, el enfoque de resolución, los resultados obtenidos y su interpretación para cada uno de los tres ejercicios del examen.

---

## Ejercicio 1 — Monte Carlo e Importance Sampling (4 puntos)

### Enunciado

Estimar $I = \int_0^1 \cos\!\left(\frac{\pi x}{2}\right) dx$ usando tres métodos: analítico, simulación directa (Monte Carlo crudo) e *importance sampling* con la densidad auxiliar $\tilde{f}(x) = \lambda(1-x^2)$.

---

### Apartado 1.1 — Determinación de $\lambda$

**Condición:** Para que $\tilde{f}$ sea densidad de probabilidad en $(0,1)$, debe integrar a 1:

$$\int_0^1 \lambda(1-x^2)\,dx = 1 \implies \lambda\left[x - \frac{x^3}{3}\right]_0^1 = 1 \implies \lambda \cdot \frac{2}{3} = 1$$

$$\boxed{\lambda = \frac{3}{2}}$$

**Por qué esta elección es buena:** La función $\tilde{f}(x) = \frac{3}{2}(1-x^2)$ tiene exactamente la misma forma cualitativa que $g(x) = \cos(\pi x/2)$: ambas valen 1 en $x=0$, se anulan en $x=1$ y son decrecientes y convexas en $(0,1)$. Cuanto más parecida sea la densidad auxiliar al integrando, menor será la varianza del estimador IS.

**Verificación numérica:** `scipy.integrate.quad` confirma $\int_0^1 \tilde{f}(x)\,dx = 1.0000000000$.

---

### Apartado 1.2 — Cálculo de $I$ por tres métodos

#### Método 1: Valor analítico exacto

$$I = \int_0^1 \cos\!\left(\frac{\pi x}{2}\right) dx = \left[\frac{2}{\pi}\sin\!\left(\frac{\pi x}{2}\right)\right]_0^1 = \frac{2}{\pi}(1 - 0) = \frac{2}{\pi}$$

$$\boxed{I = \frac{2}{\pi} \approx 0.63661977}$$

#### Método 2: Monte Carlo directo

Se interpreta la integral como esperanza bajo $\mathcal{U}(0,1)$:

$$I = \mathbb{E}_{X\sim\mathcal{U}(0,1)}[g(X)]$$

Con $N = 200$ muestras: $\hat{I}_{MC} = \frac{1}{N}\sum_{i=1}^N g(X_i)$

**Resultado obtenido:** $\hat{I}_{MC} = 0.64607647$ (error: $9.5 \times 10^{-3}$)

#### Método 3: Importance Sampling

Se reescribe la integral como esperanza bajo $\tilde{f}$:

$$I = \mathbb{E}_{X\sim\tilde{f}}\!\left[\frac{g(X)}{\tilde{f}(X)}\right] = \mathbb{E}_{X\sim\tilde{f}}\!\left[\frac{\cos(\pi X/2)}{\frac{3}{2}(1-X^2)}\right]$$

El estimador IS con $N=200$ muestras es:

$$\hat{I}_{IS} = \frac{1}{N}\sum_{i=1}^N \frac{g(X_i)}{\tilde{f}(X_i)}, \quad X_i \overset{\text{iid}}{\sim} \tilde{f}$$

**Resultado obtenido:** $\hat{I}_{IS} = 0.63715440$ (error: $5.3 \times 10^{-4}$, **18 veces menor** que el MC crudo)

#### Simulación de $\tilde{f}$: Fórmula de Cardano

Para simular de $\tilde{f}$ usamos inversión de la CDF. La CDF es:

$$F(x) = \int_0^x \frac{3}{2}(1-t^2)\,dt = \frac{3x - x^3}{2}$$

Dado $U \sim \mathcal{U}(0,1)$, resolvemos $F(x) = u$, es decir:

$$x^3 - 3x + 2u = 0$$

Esta cúbica deprimida ($p=-3$, $q=2u$) tiene discriminante $\Delta = 108(1-u^2) > 0$ para $u\in(0,1)$, lo que garantiza **tres raíces reales** (*casus irreducibilis*). Se aplica la forma trigonométrica de Cardano:

$$x_k = 2\cos\!\left(\frac{\arccos(-u)}{3} + \frac{2\pi k}{3}\right), \quad k = 0, 1, 2$$

La rama $k=1$ produce sistemáticamente la raíz en $(0,1)$ cuando $u\in(0,1)$. Verificado con un assert en el código.

---

### Apartado 1.3 — Comparación de varianzas

Para cualquier estimador Monte Carlo de la forma $\hat{I} = \frac{1}{N}\sum h(X_i)$, su varianza es:

$$\text{Var}(\hat{I}) = \frac{1}{N}\text{Var}(h(X))$$

| Método | $\text{Var}(h(X))$ | $\text{Var}(\hat{I})$ con $N=200$ |
|---|---|---|
| MC crudo | $\mathbb{E}[g^2] - I^2 = \frac{1}{2} - \frac{4}{\pi^2} \approx 0.0948$ | **0.000474** |
| Importance Sampling | $\int_0^1 g^2/\tilde{f}\,dx - I^2 \approx 0.001$ | **0.000005** |

**Reducción de varianza: 99.0%**

La varianza del IS es casi nula porque $g(x)/\tilde{f}(x) = \frac{\cos(\pi x/2)}{\frac{3}{2}(1-x^2)}$ es casi constante en $(0,1)$ (ambas funciones tienen la misma forma). Cuando el cociente es constante, la varianza es exactamente cero — este es el estimador óptimo.

**Verificación empírica:** Con 10.000 réplicas de $N=200$ muestras cada una, las varianzas empíricas son $0.000485$ (MC) y $0.000005$ (IS), coincidiendo con las teóricas.

---

## Ejercicio 2 — Transformaciones de Variables Normales (2 puntos)

### Enunciado

Sea $X\sim\mathcal{N}(0,1)$. Estudiar $Y = g(X)$ y $Z = h(X)$ con:

$$g(x) = \begin{cases}-x & x<0 \\ \alpha x & x\geq 0\end{cases}, \qquad h(x) = \begin{cases}-x & x<0 \\ \sqrt{x} & x\geq 0\end{cases}$$

---

### Apartado 2.1 — Densidad de $Y$ por simulación, análisis de $\alpha$

Se simulan $N = 100.000$ valores de $X\sim\mathcal{N}(0,1)$ y se aplica $g$ para cada $\alpha$.

| $\alpha$ | Comportamiento | Tipo de distribución | ¿Es continua? |
|---|---|---|---|
| $-1$ | $g(x) = -x$ → $Y = -X$ | $\mathcal{N}(0,1)$ idéntica | **Sí** |
| $-0.5$ | Parte positiva va a negativo | Asimétrica en $\mathbb{R}$ | **Sí** |
| $0$ | $x\geq 0$ → $Y=0$ (prob. 1/2) | **Mixta**: masa puntual en 0 + continua en $(0,\infty)$ | **No** |
| $0.5$ | Ambas ramas en $(0,\infty)$ | Unimodal asimétrica en $(0,\infty)$ | **Sí** |
| $1$ | $g(x) = |x|$ → $Y = |X|$ | **Half-normal** | **Sí** |

**$Y$ es continua si y sólo si $\alpha \neq 0$.**

- **$\alpha = 0$:** La imagen de $\{x\geq 0\}$ (prob. 1/2) se colapsa en el único punto $\{0\}$. Hay una masa puntual $P(Y=0) = 1/2$, incompatible con ser una variable continua.
- **$\alpha = 1$:** $Y = |X|$ sigue una distribución half-normal. Su densidad es $f_Y(y) = 2\phi(y)$ para $y>0$ (el doble de la normal, solo en el semieje positivo).

---

### Apartado 2.2 — Derivación analítica de $f_Y(y)$ para $\alpha > 0$

**Método: CDF → derivar**

Para $\alpha > 0$, ambas ramas de $g$ producen valores en $[0,\infty)$, por tanto $Y$ toma valores en $[0,\infty)$. Para $y > 0$:

$$P(Y \leq y) = P(\underbrace{X<0}_{\text{rama 1}}, -X\leq y) + P(\underbrace{X\geq 0}_{\text{rama 2}}, \alpha X \leq y)$$

$$= P(-y \leq X < 0) + P\!\left(0 \leq X \leq \frac{y}{\alpha}\right)$$

$$= \Phi(0) - \Phi(-y) + \Phi\!\left(\frac{y}{\alpha}\right) - \Phi(0) = \Phi\!\left(\frac{y}{\alpha}\right) - \Phi(-y)$$

Derivando respecto a $y$:

$$\boxed{f_Y(y) = \frac{1}{\alpha}\phi\!\left(\frac{y}{\alpha}\right) + \phi(y)}, \quad y > 0, \; \alpha > 0$$

**Interpretación:** La densidad es la suma de dos contribuciones gaussianas:
- $\phi(y)$: probabilidad de que la parte negativa de $X$ refleje a $y$
- $\frac{1}{\alpha}\phi(y/\alpha)$: probabilidad de que la parte positiva de $X$, escalada por $\alpha$, alcance $y$

**Verificación:** $\int_0^\infty f_Y(y)\,dy = \frac{1}{2} + \frac{1}{2} = 1$ ✓ (comprobado con `scipy.integrate.quad`)

**Casos particulares:**
- $\alpha = 1$: $f_Y(y) = \phi(y) + \phi(y) = 2\phi(y)$ → half-normal ✓
- $\alpha \to \infty$: $\frac{1}{\alpha}\phi(y/\alpha) \to 0$ → $f_Y(y) \to \phi(y)$, solo la rama negativa sobrevive

---

### Apartado 2.3 — Densidad de $Z = h(X)$: demostración y verificación

**Demostración:**

Para $z > 0$, calculamos $F_Z(z) = P(Z \leq z)$:

$$P(Z \leq z) = P(X<0,\,-X\leq z) + P(X\geq 0,\,\sqrt{X}\leq z)$$

$$= P(-z \leq X < 0) + P(0 \leq X \leq z^2)$$

$$= \Phi(0) - \Phi(-z) + \Phi(z^2) - \Phi(0) = \Phi(z^2) - \Phi(-z)$$

Derivando:

$$f_Z(z) = \frac{d}{dz}\Phi(z^2) - \frac{d}{dz}\Phi(-z) = \phi(z^2)\cdot 2z + \phi(z)$$

Como $\phi(z) = \frac{1}{\sqrt{2\pi}}e^{-z^2/2}$:

$$\boxed{f_Z(z) = \frac{1}{\sqrt{2\pi}}\left(2z\,e^{-z^4/2} + e^{-z^2/2}\right)}, \quad z > 0 \quad \blacksquare$$

**Verificación:**
- Normalización: $\int_0^\infty f_Z(z)\,dz = 1.00000000$ ✓ (con sustitución $u=z^2$ para el primer término)
- Test de Kolmogórov-Smirnov sobre $N=100.000$ muestras: $p\text{-valor} \gg 0.05$ → no se rechaza la densidad teórica ✓
- Diferencia máxima histograma–curva teórica: $< 0.01$ ✓

---

## Ejercicio 3 — Capital Económico con Cópulas (4 puntos)

### Contexto del problema

Un banco tiene exposición crediticia en 5 países. Las pérdidas se modelan como normales:

| País | $\mu_i$ | $\sigma_i$ |
|---|---|---|
| Brasil | 10 | 0.75 |
| Chile | 15 | 1.00 |
| Francia | 18 | 3.00 |
| Italia | 12 | 2.30 |
| Alemania | 19 | 4.00 |

La correlación entre países se estima a partir de datos macroeconómicos del World Bank (1991–2024): PIB, inflación, desempleo y tipo de depósito. Nivel de confianza: 95%.

---

### Paso 1 — Capital stand-alone

Con $L_i \sim \mathcal{N}(\mu_i, \sigma_i^2)$ y $z_{0.95} = 1.6449$:

$$CE_i^{\text{SA}} = \mu_i + \sigma_i \cdot 1.6449$$

| País | $CE_i^{\text{SA}}$ |
|---|---|
| Brasil | 11.2337 |
| Chile | 16.6449 |
| Francia | 22.9347 |
| Italia | 13.7834 |
| Alemania | 25.5796 |
| **TOTAL** | **92.1756** |

El capital stand-alone **asume que todos los países pierden al máximo simultáneamente** (correlación perfecta = 1). Es el escenario más conservador.

---

### Paso 2 — PCA: índice sintético macroeconómico

Se aplica PCA a las 4 series de cada país (estandarizadas, 34 observaciones anuales). El primer componente principal (PC1) maximiza la varianza explicada y actúa como barómetro del ciclo económico del país.

| País | Varianza explicada PC1 |
|---|---|
| Brasil | 59.7% |
| Chile | 59.8% |
| Francia | 55.2% |
| Italia | 37.8% |
| Alemania | 52.9% |

El PC1 de cada país resume entre el 38% y el 60% de la variabilidad conjunta de sus 4 indicadores macroeconómicos.

---

### Paso 3 — Matriz de correlaciones entre índices sintéticos

|  | Brasil | Chile | Francia | Italia | Alemania |
|---|---|---|---|---|---|
| **Brasil** | 1.000 | 0.674 | -0.205 | 0.316 | 0.024 |
| **Chile** | 0.674 | 1.000 | -0.241 | 0.562 | 0.031 |
| **Francia** | -0.205 | -0.241 | 1.000 | 0.284 | 0.544 |
| **Italia** | 0.316 | 0.562 | 0.284 | 1.000 | 0.112 |
| **Alemania** | 0.024 | 0.031 | 0.544 | 0.112 | 1.000 |

**Observaciones económicas:**
- Brasil y Chile muestran alta correlación (0.674): economías latinoamericanas con ciclos similares.
- Francia y Alemania están moderadamente correlacionadas (0.544): integración en la zona euro.
- Los países latinoamericanos tienen correlación negativa con Francia (−0.2/−0.24): ciclos divergentes.
- Alemania y los países latinoamericanos prácticamente no correlacionan (~0.03).

La matriz es semidefinida positiva (todos los valores propios ≥ 0) → válida para simulación.

---

### Paso 4 — Capital diversificado: Cópula Gaussiana

**Fundamento:** La cópula Gaussiana modela la dependencia entre los cuantiles uniformes de las pérdidas mediante una normal multivariante.

**Algoritmo (Cholesky + simulación):**

1. Factorizar: $\Sigma = L L^\top$ (descomposición de Cholesky de la matriz de correlaciones)
2. Simular $\mathbf{z} \sim \mathcal{N}_5(\mathbf{0}, I_5)$ independientes
3. Construir $\mathbf{Z} = L\mathbf{z} \sim \mathcal{N}_5(\mathbf{0}, \Sigma)$
4. Dado que las marginales son normales: $L_i = \mu_i + \sigma_i Z_i$
5. $L_{\text{total}} = \sum_i L_i$
6. $CE_{\text{G}} = \text{percentil}_{95\%}(L_{\text{total}})$ sobre $N = 500.000$ simulaciones

**Resultado:** $CE_{\text{Gaussiana}} = \mathbf{86.03}$

---

### Paso 5 — Capital diversificado: Cópula t-Student

**Motivación:** La cópula Gaussiana no captura *tail dependence* (dependencia en las colas). En crisis financieras, los países tienden a correlacionarse más de lo que predice la normal. La cópula t-Student, con colas más pesadas, modela este fenómeno.

**Algoritmo:**

1. $\mathbf{Z} \sim \mathcal{N}_5(\mathbf{0}, \Sigma)$ correlacionada (como antes)
2. $\chi^2 \sim \chi^2_\nu$ independiente (con $\nu = 5$ grados de libertad)
3. $\mathbf{T} = \mathbf{Z} / \sqrt{\chi^2/\nu}$ → t multivariante con $\nu$ g.l.
4. $U_i = F_{t_\nu}(T_i)$ → uniformes con dependencia t
5. $L_i = \mu_i + \sigma_i \Phi^{-1}(U_i)$ → marginales normales con cópula t
6. $CE_{\text{t}} = \text{percentil}_{95\%}(L_{\text{total}})$

**Resultado:** $CE_{\text{t-Student}} = \mathbf{85.93}$

> **Nota:** En este caso concreto, el capital de la cópula t es ligeramente inferior al Gaussiano. Esto puede parecer contraintuitivo (la t tiene colas más pesadas), pero depende de la estructura de correlación específica estimada. Con correlaciones bajas o negativas entre algunos pares (como los países latinoamericanos con los europeos), la mayor dependencia en las colas de la cópula t puede no aumentar el percentil 95% del total.

---

### Paso 6 — Comparación y discusión

| Método | Capital (p95) | Beneficio div. | % ahorro |
|---|---|---|---|
| Stand-alone (sin div.) | 92.18 | — | — |
| Cópula Gaussiana | 86.03 | 6.15 | 6.7% |
| Cópula t-Student (ν=5) | 85.93 | 6.25 | 6.8% |

**Impacto de la diversificación:**
El beneficio de ~6.7% proviene de que los países no están todos perfectamente correlacionados. Las pérdidas máximas individuales no ocurren todas al mismo tiempo, permitiendo al banco liberar capital.

**Diferencia entre cópulas:**
La diferencia entre Gaussiana y t es pequeña en este caso (≈0.10 u.m.) debido a las correlaciones relativamente bajas entre los países (especialmente los latinoamericanos con los europeos). Si las correlaciones fuesen más altas, la cópula t produciría un capital notablemente mayor al capturar la mayor co-dependencia en las colas.

**Conclusión regulatoria:**
En general, la cópula t-Student es más prudente para modelar escenarios de estrés sistémico. El impacto depende de la estructura de correlación y de los grados de libertad elegidos (menor ν = colas más pesadas = mayor capital).

---

*Fin de la memoria de resolución.*
