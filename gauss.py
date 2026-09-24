"""
Resolución de sistemas de ecuaciones lineales por el método de Gauss.

El sistema tiene m ecuaciones y n variables, con m <= n, y se introduce como
una matriz ampliada rectangular de m filas y n + 1 columnas:

    a11 x1 + a12 x2 + ... + a1n xn = b1          [ a11 a12 ... a1n | b1 ]
    a21 x1 + a22 x2 + ... + a2n xn = b2   --->   [ a21 a22 ... a2n | b2 ]
    ...                                          [ ...             | .. ]
    am1 x1 + am2 x2 + ... + amn xn = bm          [ am1 am2 ... amn | bm ]

El método de Gauss consta de dos fases:
  1. Eliminación hacia adelante: con operaciones elementales de fila se lleva
     la matriz a forma escalonada.
  2. Sustitución hacia atrás: se despejan las variables empezando por la
     última ecuación no nula.

Se usa aritmética exacta con fracciones (fractions.Fraction) para evitar
errores de redondeo. Según el rango de la matriz de coeficientes (r) y el de
la matriz ampliada (r'), el sistema se clasifica (teorema de Rouché-Frobenius):
  - r != r'           -> incompatible (sin solución)
  - r == r' == n      -> compatible determinado (solución única)
  - r == r' <  n      -> compatible indeterminado (infinitas soluciones,
                         expresadas en función de n - r parámetros)
"""

from fractions import Fraction

INCOMPATIBLE = "incompatible"
DETERMINADO = "compatible determinado"
INDETERMINADO = "compatible indeterminado"


def a_fraccion(valor):
    """Convierte un número o texto ("3", "-2.5", "1/3") a Fraction."""
    if isinstance(valor, Fraction):
        return valor
    if isinstance(valor, float):
        return Fraction(str(valor))
    return Fraction(valor)


def validar_matriz(matriz):
    """Comprueba que la matriz ampliada sea válida y devuelve (m, n)."""
    m = len(matriz)
    if m == 0:
        raise ValueError("La matriz no puede estar vacía.")
    columnas = len(matriz[0])
    if columnas < 2:
        raise ValueError("La matriz debe tener al menos 2 columnas "
                         "(una variable y el término independiente).")
    for i, fila in enumerate(matriz, start=1):
        if len(fila) != columnas:
            raise ValueError(f"La fila {i} tiene {len(fila)} elementos; "
                             f"se esperaban {columnas}.")
    n = columnas - 1
    if m > n:
        raise ValueError(f"El número de ecuaciones (m = {m}) debe ser menor "
                         f"o igual que el número de variables (n = {n}).")
    return m, n


def eliminacion_gauss(matriz, pasos=None):
    """
    Fase 1: lleva la matriz ampliada a forma escalonada por filas.

    Usa pivoteo parcial (se elige como pivote el elemento de mayor valor
    absoluto de la columna) e intercambio de filas cuando es necesario.
    Devuelve (matriz_escalonada, columnas_pivote). Si `pasos` es una lista,
    se le añade una descripción de cada operación realizada.
    """
    m, n = validar_matriz(matriz)
    a = [[a_fraccion(x) for x in fila] for fila in matriz]
    columnas_pivote = []
    fila = 0

    for col in range(n):
        if fila >= m:
            break
        # Pivoteo parcial: fila con el mayor |a[i][col]| desde `fila`.
        p = max(range(fila, m), key=lambda i: abs(a[i][col]))
        if a[p][col] == 0:
            continue  # columna sin pivote: corresponde a una variable libre
        if p != fila:
            a[fila], a[p] = a[p], a[fila]
            if pasos is not None:
                pasos.append((f"F{fila + 1} <-> F{p + 1}", copiar(a)))
        for i in range(fila + 1, m):
            if a[i][col] != 0:
                factor = a[i][col] / a[fila][col]
                for j in range(col, n + 1):
                    a[i][j] -= factor * a[fila][j]
                if pasos is not None:
                    pasos.append((f"F{i + 1} <- F{i + 1} - ({factor})·F{fila + 1}",
                                  copiar(a)))
        columnas_pivote.append(col)
        fila += 1

    return a, columnas_pivote


def sustitucion_hacia_atras(a, columnas_pivote, n):
    """
    Fase 2: resuelve el sistema escalonado por sustitución hacia atrás.

    Cada variable se representa como una expresión afín de los parámetros
    libres: un diccionario {índice_variable_libre: coeficiente, None: constante}.
    """
    libres = [j for j in range(n) if j not in columnas_pivote]
    solucion = [None] * n
    for j in libres:
        solucion[j] = {j: Fraction(1), None: Fraction(0)}

    for i in reversed(range(len(columnas_pivote))):
        c = columnas_pivote[i]
        pivote = a[i][c]
        expr = {None: a[i][n] / pivote}
        for j in range(c + 1, n):
            if a[i][j] == 0:
                continue
            k = -a[i][j] / pivote
            for clave, valor in solucion[j].items():
                expr[clave] = expr.get(clave, Fraction(0)) + k * valor
        solucion[c] = {k: v for k, v in expr.items() if v != 0 or k is None}

    return solucion, libres


def resolver_gauss(matriz, pasos=None):
    """
    Resuelve el sistema dado por su matriz ampliada m x (n + 1), con m <= n.

    Devuelve un diccionario con:
      - "tipo":      INCOMPATIBLE, DETERMINADO o INDETERMINADO
      - "rango":     rango de la matriz de coeficientes
      - "escalonada": matriz escalonada obtenida
      - "solucion":  None si es incompatible; en otro caso, una lista con la
                     expresión de cada variable (ver sustitucion_hacia_atras)
      - "libres":    índices de las variables libres (parámetros)
    """
    m, n = validar_matriz(matriz)
    a, columnas_pivote = eliminacion_gauss(matriz, pasos)
    rango = len(columnas_pivote)

    resultado = {"rango": rango, "escalonada": a, "solucion": None, "libres": []}

    # Una fila 0 0 ... 0 | b con b != 0 indica un sistema incompatible.
    for i in range(rango, m):
        if a[i][n] != 0:
            resultado["tipo"] = INCOMPATIBLE
            return resultado

    solucion, libres = sustitucion_hacia_atras(a, columnas_pivote, n)
    resultado["solucion"] = solucion
    resultado["libres"] = libres
    resultado["tipo"] = DETERMINADO if rango == n else INDETERMINADO
    return resultado


# ---------------------------------------------------------------------------
# Presentación de resultados
# ---------------------------------------------------------------------------

def copiar(a):
    return [fila[:] for fila in a]


def formato_numero(x):
    x = a_fraccion(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def formato_matriz(a):
    """Devuelve la matriz ampliada como texto, con la barra separadora."""
    textos = [[formato_numero(x) for x in fila] for fila in a]
    ancho = max(len(t) for fila in textos for t in fila)
    lineas = []
    for fila in textos:
        izquierda = " ".join(t.rjust(ancho) for t in fila[:-1])
        lineas.append(f"[ {izquierda} | {fila[-1].rjust(ancho)} ]")
    return "\n".join(lineas)


def formato_expresion(expr, parametros):
    """Convierte {j: coef, None: cte} en texto, p. ej. '2 - 3·t1 + t2'."""
    partes = []
    constante = expr.get(None, Fraction(0))
    for j in sorted(k for k in expr if k is not None):
        coef = expr[j]
        if coef == 0:
            continue
        signo = "-" if coef < 0 else "+"
        valor = abs(coef)
        termino = parametros[j] if valor == 1 else f"{formato_numero(valor)}·{parametros[j]}"
        partes.append((signo, termino))

    if constante != 0 or not partes:
        texto = formato_numero(constante)
    else:
        signo, termino = partes.pop(0)
        texto = termino if signo == "+" else f"-{termino}"
    for signo, termino in partes:
        texto += f" {signo} {termino}"
    return texto


def formato_solucion(resultado):
    """Texto legible con la clasificación y la solución del sistema."""
    tipo = resultado["tipo"]
    if tipo == INCOMPATIBLE:
        return "Sistema INCOMPATIBLE: no tiene solución."

    solucion = resultado["solucion"]
    libres = resultado["libres"]
    parametros = {j: f"t{k}" for k, j in enumerate(libres, start=1)}
    lineas = []
    if tipo == DETERMINADO:
        lineas.append("Sistema COMPATIBLE DETERMINADO: solución única.")
    else:
        lineas.append(f"Sistema COMPATIBLE INDETERMINADO: infinitas soluciones "
                      f"({len(libres)} parámetro(s) libre(s)).")
    for j, expr in enumerate(solucion):
        if j in parametros:
            lineas.append(f"  x{j + 1} = {parametros[j]}   (parámetro libre)")
        else:
            lineas.append(f"  x{j + 1} = {formato_expresion(expr, parametros)}")
    return "\n".join(lineas)


# ---------------------------------------------------------------------------
# Entrada interactiva
# ---------------------------------------------------------------------------

def leer_entero(mensaje, minimo=1):
    while True:
        try:
            valor = int(input(mensaje))
            if valor >= minimo:
                return valor
            print(f"  Debe ser un entero mayor o igual que {minimo}.")
        except ValueError:
            print("  Entrada no válida, introduce un número entero.")


def leer_matriz():
    """Pide al usuario m, n y la matriz ampliada fila por fila."""
    n = leer_entero("Número de variables (n): ")
    while True:
        m = leer_entero("Número de ecuaciones (m, con m <= n): ")
        if m <= n:
            break
        print(f"  El número de ecuaciones debe ser <= {n}.")

    print(f"\nIntroduce la matriz ampliada de {m} filas y {n + 1} columnas.")
    print("Escribe cada fila con los valores separados por espacios "
          "(se admiten enteros, decimales y fracciones como 1/3).")
    matriz = []
    for i in range(1, m + 1):
        while True:
            texto = input(f"  Fila {i}: ").replace(",", " ").split()
            if len(texto) != n + 1:
                print(f"  Se esperaban {n + 1} valores y se introdujeron {len(texto)}.")
                continue
            try:
                matriz.append([a_fraccion(x) for x in texto])
                break
            except (ValueError, ZeroDivisionError):
                print("  Algún valor no es un número válido.")
    return matriz


def main():
    print("=" * 60)
    print(" Resolución de sistemas lineales por el método de Gauss")
    print("=" * 60)
    matriz = leer_matriz()

    print("\nMatriz ampliada inicial:")
    print(formato_matriz(matriz))

    pasos = []
    resultado = resolver_gauss(matriz, pasos)

    if pasos:
        print("\nEliminación hacia adelante:")
        for descripcion, estado in pasos:
            print(f"\n{descripcion}")
            print(formato_matriz(estado))

    print("\nMatriz escalonada:")
    print(formato_matriz(resultado["escalonada"]))
    print(f"\nRango de la matriz de coeficientes: {resultado['rango']}")
    print()
    print(formato_solucion(resultado))


if __name__ == "__main__":
    main()
