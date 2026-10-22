# Matematicas-GAUSS

## App para el teléfono: Gauss |A|

La carpeta [`app/`](app/) contiene una app web instalable (PWA) que funciona en el iPhone:

- **Determinantes** de matrices de 1×1 a 8×8 por eliminación de Gauss, con comprobación por
  `ad − bc` (2×2) y por la regla de Sarrus (3×3).
- **Sistemas de ecuaciones** con m ecuaciones y n variables (m ≤ n): solución única, infinitas
  soluciones con parámetros o sistema incompatible (Rouché–Frobenius).
- **Dashboard** con el resultado arriba. La matriz se **genera automáticamente** (rango de valores y
  fracciones opcionales) o se **escribe a mano**, y el resultado se recalcula al escribir.
- Paso a paso con cada operación de fila y el pivote resaltado, fracciones exactas, pivoteo parcial
  opcional, historial y funcionamiento sin conexión.

### Instalarla en el iPhone

1. Publica la app: en GitHub ve a **Settings → Pages** y en **Source** elige **GitHub Actions**.
   Cada cambio en `main` la vuelve a publicar en `https://<usuario>.github.io/Matematicas-GAUSS/`.
2. Abre esa dirección en **Safari**, toca **Compartir** y luego **Añadir a pantalla de inicio**.

Para probarla en tu computadora: `cd app && python3 -m http.server 8000` y abre http://localhost:8000.

Pruebas de la app: `node --test tests/gauss.test.js`

## Programa de consola (Python)

Programa en Python que resuelve sistemas de ecuaciones lineales por el **método de Gauss**.

El sistema tiene **n variables** y **m ecuaciones**, con **m ≤ n**, y se introduce como una
matriz ampliada rectangular de **m filas** y **n + 1 columnas** (coeficientes y término independiente).

### Uso

```bash
python3 gauss.py
```

El programa pide `n`, `m` y cada fila de la matriz (se admiten enteros, decimales y fracciones como `1/3`).
Muestra cada operación de fila de la eliminación, la matriz escalonada, el rango y la solución:

- **Compatible determinado**: solución única.
- **Compatible indeterminado**: infinitas soluciones en función de parámetros `t1, t2, ...`.
- **Incompatible**: sin solución.

También puede usarse como módulo:

```python
from gauss import resolver_gauss, formato_solucion

resultado = resolver_gauss([[2, 1, -1, 8],
                            [-3, -1, 2, -11],
                            [-2, 1, 2, -3]])
print(formato_solucion(resultado))   # x1 = 2, x2 = 3, x3 = -1
```

Los cálculos se hacen con fracciones exactas, sin errores de redondeo.

### Pruebas

```bash
python3 -m unittest -v test_gauss
```
