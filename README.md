# Matematicas-GAUSS

Programa en Python que resuelve sistemas de ecuaciones lineales por el **método de Gauss**.

El sistema tiene **n variables** y **m ecuaciones**, con **m ≤ n**, y se introduce como una
matriz ampliada rectangular de **m filas** y **n + 1 columnas** (coeficientes y término independiente).

## Uso

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

## Pruebas

```bash
python3 -m unittest -v test_gauss
```
