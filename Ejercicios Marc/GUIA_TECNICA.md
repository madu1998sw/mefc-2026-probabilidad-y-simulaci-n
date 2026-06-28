# Guía Técnica del Proyecto
## Máster Executive en Finanzas Cuantitativas 2026 — AFI Global Education

---

> Este documento explica **cómo está organizado el proyecto**, qué decisiones técnicas se tomaron y por qué, y cómo reproducir todos los resultados de principio a fin. Es el complemento técnico de la `MEMORIA_RESOLUCION.md`.

---

## 1. Estructura del repositorio

```
MEFC-2026-PROBABILIDAD-Y-SIMULACION/
│
├── 0-enunciado/
│   └── MEFC_2026_examen_probabilidad_y_simulacion.pdf
│
├── 1-monte-carlo-importance-sampling/
│   ├── ejercicio1_monte_carlo_importance_sampling.ipynb  ← ENTREGABLE
│   └── resultados/
│       ├── grafico_01_integrando_vs_ftilde.png
│       ├── grafico_02_estimadores.png
│       └── grafico_03_varianzas_estimadores.png
│
├── 2-transformaciones-normales/
│   ├── ejercicio2_transformaciones_normales.ipynb        ← ENTREGABLE
│   └── resultados/
│       ├── grafico_04_fY_por_alpha.png
│       ├── grafico_05_fY_analitica_vs_sim.png
│       └── grafico_06_fZ_verificacion.png
│
├── 3-capital-economico-copulas/
│   ├── ejercicio3_capital_economico_copulas.ipynb        ← ENTREGABLE
│   └── resultados/
│       ├── grafico_07_capital_standalone.png
│       ├── grafico_08_indices_pca.png
│       ├── grafico_09_heatmap_correlaciones.png
│       └── grafico_10_comparacion_capitales.png
│
├── series_macro__1_.xlsx                                 ← DATOS (World Bank)
├── MEMORIA_RESOLUCION.md                                 ← este archivo
├── GUIA_TECNICA.md                                       ← guía técnica
└── README.md
```

**Regla fundamental:** cada ejercicio vive en su propia carpeta, con su notebook y su subcarpeta `resultados/`. Los notebooks son el entregable principal; los `.png` son evidencia visual generada al ejecutarlos.

---

## 2. Enfoque general de resolución

### Filosofía

El trabajo sigue un enfoque de **verificación cruzada**: cada resultado se obtiene analíticamente (cuando es posible) y se confirma mediante simulación. Si hay discrepancia, el error está en el código, no en las matemáticas.

Tres niveles de verificación por resultado:
1. **Desarrollo analítico** en celdas Markdown con LaTeX
2. **Implementación en Python** con comentarios paso a paso
3. **Confirmación numérica/gráfica**: `scipy.integrate.quad`, test KS, diferencia histograma–curva teórica

### Por qué notebooks Jupyter y no scripts `.py`

Los notebooks permiten intercalar matemáticas (LaTeX), código Python y gráficos en un único documento lineal. Esto es ideal para una memoria académica: el profesor puede ver el razonamiento, el código y el resultado en el mismo artefacto, en el orden en que se desarrolla.

---

## 3. Convenciones de código comunes a todos los notebooks

### Setup centralizado (primera celda de código)

Cada notebook comienza con una única celda que contiene:

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate
# ... más imports específicos del ejercicio

plt.rcParams.update({...})   # estilo coherente
SEED = 42
rng  = np.random.default_rng(SEED)  # generador reproducible
N    = 200   # constantes del enunciado
os.makedirs("resultados", exist_ok=True)
```

**Por qué una sola celda de setup:** Si los imports, la semilla y las constantes están dispersos por el notebook, un cambio en un parámetro requiere editar múltiples celdas. Centralizado, basta cambiar una línea.

### Generador de números aleatorios moderno

Se usa `np.random.default_rng(SEED)` (API moderna de NumPy, disponible desde v1.17) en lugar del antiguo `np.random.seed()`. Ventajas:
- El estado del generador es explícito (no global)
- Permite pasar generadores distintos a funciones independientes sin interferencia
- Mayor calidad estadística (algoritmo PCG64)

### Aserciones defensivas

```python
assert abs(integral - 1.0) < 1e-10, "La densidad no integra a 1"
assert np.all(eigenvalues >= -1e-8), "Matriz no semidefinida positiva"
```

Las aserciones fallan en voz alta si hay un error de implementación. Son especialmente importantes en:
- Verificación de densidades (deben integrar a 1)
- Selección de la raíz correcta en la cúbica de Cardano
- Validez de la matriz de correlaciones para Cholesky

### Figuras: siempre teoría sobre simulación

Cada figura que compara simulación con teoría sigue el mismo patrón:
```python
ax.hist(muestras, density=True, ...)   # histograma simulado
ax.plot(x, f_teorica(x), ...)          # curva analítica encima
```

Esto permite verificar visualmente en un golpe de vista si la simulación es correcta.

---

## 4. Ejercicio 1: decisiones técnicas

### Selección de la raíz de Cardano

La cúbica $x^3 - 3x + 2u = 0$ tiene tres raíces reales para $u \in (0,1)$. La solución trigonométrica produce:

```python
theta = np.arccos(-u) / 3
roots = np.stack([
    2 * np.cos(theta),
    2 * np.cos(theta + 2*np.pi/3),   # ← esta es la rama en (0,1)
    2 * np.cos(theta + 4*np.pi/3),
], axis=1)
mask = (roots > 0) & (roots < 1)
samples = roots[mask].reshape(n)
```

Se verifica con `assert` que exactamente una raíz por muestra cae en $(0,1)$.

### Cálculo de varianzas

- **MC crudo:** $\mathbb{E}[g^2] = \int_0^1 \cos^2(\pi x/2)\,dx = 1/2$ (analítico exacto por identidad trigonométrica)
- **IS:** $\mathbb{E}[w^2] = \int_0^1 g(x)^2/\tilde{f}(x)\,dx$ calculado con `scipy.integrate.quad` (sin forma cerrada)
- El denominador de IS se evalúa en $[0, 1-10^{-10}]$ para evitar la singularidad en $x=1$ donde $\tilde{f}(1) = 0$

---

## 5. Ejercicio 2: decisiones técnicas

### $N = 100.000$ muestras

El enunciado no especifica $N$ para las simulaciones del Ejercicio 2. Se usaron 100.000 para obtener histogramas suaves y un test KS fiable. Más muestras no añaden información útil en este contexto.

### Caso $\alpha = 0$: variable mixta

Al simular con $\alpha = 0$, exactamente $\approx 50.000$ muestras colapsan a $Y = 0$. El histograma ordinario no muestra correctamente la masa puntual (la barra en 0 es infinitamente alta en densidad). En el código se usa un histograma con `range=(lo, hi)` excluyendo la masa puntual para visualizar la parte continua, y se comenta en el Markdown.

### Test de Kolmogórov-Smirnov para $f_Z$

```python
ks_stat, ks_pval = stats.kstest(Z, F_Z_analitica)
```

El KS compara la CDF empírica de las muestras con la CDF teórica $F_Z(z) = \Phi(z^2) - \Phi(-z)$. Un p-valor alto confirma que la simulación sigue la densidad derivada analíticamente.

---

## 6. Ejercicio 3: decisiones técnicas

### Lectura del Excel (`series_macro.xlsx`)

```python
df_raw = pd.read_excel(EXCEL_PATH, sheet_name='Data')
year_cols = [c for c in df_raw.columns if str(c)[:4].isdigit()]
```

Las columnas de años tienen el formato `"1991 [YR1991]"`, por eso se filtran por los primeros 4 dígitos. Hay 34 años (1991–2024), 5 países × 4 series = 20 filas en total.

### PCA con `n_components=1`

Se especifica `n_components=1` (no 4) porque solo necesitamos el primer componente. Usar 4 con datos 34×4 causaría un error en algunas versiones de sklearn cuando el número de componentes iguala el mínimo de dimensiones.

```python
pca = PCA(n_components=1, random_state=SEED)
pc1 = pca.fit_transform(X_std)[:, 0]
```

### Cholesky para la normal multivariante

La simulación de $\mathcal{N}_5(\mathbf{0}, \Sigma)$ se hace mediante:

```python
L = np.linalg.cholesky(corr_matrix)   # Σ = L·Lᵀ
Z = rng.standard_normal((N_SIM, 5)) @ L.T
```

Es más eficiente que `np.random.multivariate_normal` para muestras grandes porque evita recomputar la factorización en cada llamada.

### Semillas separadas para Gaussiana y t

```python
rng_gauss = np.random.default_rng(SEED)
rng_t     = np.random.default_rng(SEED + 10)
```

Usar semillas distintas garantiza que ambas simulaciones son independientes entre sí, evitando correlaciones artificiales en la comparación de resultados.

### Número de simulaciones: $N = 500.000$

El percentil 95% sobre una muestra de $N$ observaciones tiene un error estándar aproximado de:

$$\text{SE}(CE) \approx \frac{\sigma_L}{f_L(CE)\sqrt{N}}$$

Con $N = 500.000$, este error es del orden de $0.01$ u.m., suficiente para comparar los tres métodos. Con $N = 10.000$ la estimación del percentil sería demasiado ruidosa.

---

## 7. Cómo ejecutar los notebooks

### Requisitos

```bash
pip install numpy pandas matplotlib scipy scikit-learn openpyxl jupyterlab
```

### Ejecución completa (recomendado)

Desde la raíz del proyecto, ejecutar cada notebook:

```bash
# Ejercicio 1
jupyter nbconvert --to notebook --execute --inplace \
  1-monte-carlo-importance-sampling/ejercicio1_monte_carlo_importance_sampling.ipynb

# Ejercicio 2
jupyter nbconvert --to notebook --execute --inplace \
  2-transformaciones-normales/ejercicio2_transformaciones_normales.ipynb

# Ejercicio 3
jupyter nbconvert --to notebook --execute --inplace \
  3-capital-economico-copulas/ejercicio3_capital_economico_copulas.ipynb
```

### En VSCode

1. Abrir la carpeta raíz del proyecto
2. Abrir el notebook `.ipynb` deseado
3. Seleccionar el kernel de Python del entorno virtual
4. Pulsar **Run All** (Ctrl+Shift+P → "Run All Cells")
5. Las figuras se guardan automáticamente en `resultados/`

### Nota sobre la ruta del Excel (Ejercicio 3)

La ruta en el notebook es relativa:
```python
EXCEL_PATH = '../../../series_macro__1_.xlsx'
```
Esto asume que el Excel está en la raíz del proyecto y el notebook está 3 niveles por debajo. **Ajustar si la estructura es diferente.**

---

## 8. Resultados exactos por ejercicio

### Ejercicio 1

| Métrica | Valor |
|---|---|
| $\lambda$ | 1.5 |
| $I$ analítico | 0.63661977 |
| $\hat{I}_{MC}$ (seed=42, N=200) | 0.64607647 |
| $\hat{I}_{IS}$ (seed=42, N=200) | 0.63715440 |
| $\text{Var}(\hat{I}_{MC})$ | 0.000474 |
| $\text{Var}(\hat{I}_{IS})$ | 0.000005 |
| Reducción de varianza IS vs MC | **99.0%** |

### Ejercicio 2

| Métrica | Valor |
|---|---|
| $Y$ continua | $\alpha \neq 0$ |
| $\alpha=0$: tipo | Mixta (masa puntual en 0) |
| $\alpha=1$: tipo | Half-normal |
| $\int_0^\infty f_Z(z)\,dz$ | 1.00000000 |
| KS p-valor ($f_Z$) | $\gg 0.05$ (no se rechaza) |

### Ejercicio 3

| Métrica | Valor |
|---|---|
| Capital stand-alone total | 92.1756 |
| Capital diversificado (Gaussiana) | 86.03 |
| Capital diversificado (t, ν=5) | 85.93 |
| Beneficio div. (Gaussiana) | 6.15 (6.7%) |
| Beneficio div. (t) | 6.25 (6.8%) |

---

## 9. Checklist de entrega

- [x] Ejercicio 1: notebook completo, 3 gráficos, todos los outputs ejecutados
- [x] Ejercicio 2: notebook completo, 3 gráficos, test KS, todos los outputs ejecutados
- [x] Ejercicio 3: notebook completo, 4 gráficos, todos los outputs ejecutados
- [x] `MEMORIA_RESOLUCION.md`: desarrollo teórico completo A–Z
- [x] `GUIA_TECNICA.md`: este documento
- [x] Todos los resultados numéricos de las memorias coinciden con los outputs ejecutados
- [x] Semilla fija `SEED = 42` en todos los notebooks → reproducibilidad garantizada
- [x] `series_macro__1_.xlsx` incluido en el repositorio
