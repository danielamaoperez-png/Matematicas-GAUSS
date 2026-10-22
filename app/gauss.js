/*
 * Lógica matemática de la app: fracciones exactas, determinantes por el
 * método de Gauss y resolución de sistemas lineales m x n (m <= n).
 * Es la misma lógica que gauss.py, traducida a JavaScript.
 */
(function (global) {
  "use strict";

  // ------------------------------------------------------------------
  // Fracciones exactas (numerador y denominador BigInt)
  // ------------------------------------------------------------------
  function mcd(a, b) {
    if (a < 0n) a = -a;
    if (b < 0n) b = -b;
    while (b) [a, b] = [b, a % b];
    return a;
  }

  class Q {
    constructor(n, d = 1n) {
      if (d === 0n) throw new Error("División entre cero.");
      if (d < 0n) { n = -n; d = -d; }
      const g = mcd(n, d) || 1n;
      this.n = n / g;
      this.d = d / g;
    }
    add(o) { return new Q(this.n * o.d + o.n * this.d, this.d * o.d); }
    sub(o) { return new Q(this.n * o.d - o.n * this.d, this.d * o.d); }
    mul(o) { return new Q(this.n * o.n, this.d * o.d); }
    div(o) {
      if (o.n === 0n) throw new Error("División entre cero.");
      return new Q(this.n * o.d, this.d * o.n);
    }
    neg() { return new Q(-this.n, this.d); }
    abs() { return new Q(this.n < 0n ? -this.n : this.n, this.d); }
    esCero() { return this.n === 0n; }
    esUno() { return this.n === 1n && this.d === 1n; }
    signo() { return this.n < 0n ? -1 : this.n > 0n ? 1 : 0; }
    cmp(o) {
      const l = this.n * o.d, r = o.n * this.d;
      return l < r ? -1 : l > r ? 1 : 0;
    }
    eq(o) { return this.n === o.n && this.d === o.d; }
    toString() { return this.d === 1n ? `${this.n}` : `${this.n}/${this.d}`; }
    aNumero() { return Number(this.n) / Number(this.d); }
  }

  const CERO = new Q(0n);
  const UNO = new Q(1n);

  function decimalAQ(texto) {
    const m = /^([+-]?)(\d*)(?:\.(\d*))?$/.exec(texto);
    if (!m || (m[2] === "" && (m[3] === undefined || m[3] === ""))) return null;
    const enteros = m[2] || "0";
    const dec = m[3] || "";
    const q = new Q(BigInt(enteros + dec), 10n ** BigInt(dec.length));
    return m[1] === "-" ? q.neg() : q;
  }

  /**
   * Convierte un texto ("3", "-2.5", "2,5", "1/3", "−4") a fracción.
   * Una celda vacía vale 0. Devuelve null si el texto no es un número.
   */
  function parsear(texto) {
    if (texto instanceof Q) return texto;
    if (typeof texto === "number") texto = String(texto);
    const t = String(texto).trim().replace(/\s+/g, "").replace(/,/g, ".").replace(/[−–]/g, "-");
    if (t === "") return CERO;
    const partes = t.split("/");
    if (partes.length === 1) return decimalAQ(partes[0]);
    if (partes.length !== 2) return null;
    const a = decimalAQ(partes[0]), b = decimalAQ(partes[1]);
    if (!a || !b || b.esCero()) return null;
    return a.div(b);
  }

  /** Convierte una matriz de textos/números a fracciones; lanza error con la celda mala. */
  function aMatrizQ(filas) {
    return filas.map((fila, i) => fila.map((v, j) => {
      const q = parsear(v);
      if (!q) {
        const e = new Error(`El valor «${v}» de la fila ${i + 1}, columna ${j + 1} no es un número válido.`);
        e.celda = [i, j];
        throw e;
      }
      return q;
    }));
  }

  const copiar = (a) => a.map((f) => f.slice());

  function elegirPivote(a, desde, m, col, pivoteoParcial) {
    let p = -1;
    for (let i = desde; i < m; i++) {
      if (a[i][col].esCero()) continue;
      if (p < 0) { p = i; if (!pivoteoParcial) break; continue; }
      if (a[i][col].abs().cmp(a[p][col].abs()) > 0) p = i;
    }
    return p;
  }

  // ------------------------------------------------------------------
  // Determinante por eliminación de Gauss
  // ------------------------------------------------------------------
  /**
   * Calcula det(A) llevando A a forma triangular superior.
   *   - Intercambiar dos filas cambia el signo del determinante.
   *   - Sumar a una fila un múltiplo de otra no lo cambia.
   *   - det de una triangular = producto de la diagonal.
   */
  function determinante(entrada, opciones = {}) {
    const pivoteoParcial = opciones.pivoteoParcial !== false;
    const original = aMatrizQ(entrada);
    const n = original.length;
    if (n === 0 || original.some((f) => f.length !== n)) {
      throw new Error("La matriz debe ser cuadrada (mismo número de filas y columnas).");
    }
    const a = copiar(original);
    const pasos = [];
    let intercambios = 0;

    for (let c = 0; c < n; c++) {
      const p = elegirPivote(a, c, n, c, pivoteoParcial);
      if (p < 0) {
        pasos.push({
          op: `Columna ${c + 1} sin pivote`,
          nota: `Todos los elementos de la columna ${c + 1} desde la fila ${c + 1} son 0, así que la matriz es singular y det(A) = 0.`,
          m: copiar(a), columna: c,
        });
        return { det: CERO, singular: true, intercambios, diagonal: null, triangular: a, pasos, original, n };
      }
      if (p !== c) {
        [a[c], a[p]] = [a[p], a[c]];
        intercambios++;
        pasos.push({
          op: `F${c + 1} ↔ F${p + 1}`,
          nota: "Intercambiar dos filas cambia el signo del determinante.",
          m: copiar(a), piv: [c, c], filas: [c, p], tipo: "swap",
        });
      }
      for (let i = c + 1; i < n; i++) {
        if (a[i][c].esCero()) continue;
        const f = a[i][c].div(a[c][c]);
        for (let j = c; j < n; j++) a[i][j] = a[i][j].sub(f.mul(a[c][j]));
        pasos.push({
          op: operacionFila(i, f, c),
          nota: "Sumar a una fila un múltiplo de otra no cambia el determinante.",
          m: copiar(a), piv: [c, c], filas: [i],
        });
      }
    }

    const diagonal = a.map((f, i) => f[i]);
    let det = diagonal.reduce((acc, x) => acc.mul(x), UNO);
    if (intercambios % 2) det = det.neg();
    return { det, singular: det.esCero(), intercambios, diagonal, triangular: a, pasos, original, n };
  }

  /** Comprobación por fórmula directa para 1x1, 2x2 (ad − bc) y 3x3 (regla de Sarrus). */
  function comprobacionDirecta(A) {
    const n = A.length;
    if (n === 1) return { metodo: "Matriz 1×1", valor: A[0][0], expresion: `det(A) = ${A[0][0]}` };
    if (n === 2) {
      const [[a, b], [c, d]] = A;
      return {
        metodo: "Fórmula 2×2: ad − bc",
        valor: a.mul(d).sub(b.mul(c)),
        expresion: `(${a})(${d}) − (${b})(${c})`,
      };
    }
    if (n === 3) {
      const [[a, b, c], [d, e, f], [g, h, i]] = A;
      const pos = a.mul(e).mul(i).add(b.mul(f).mul(g)).add(c.mul(d).mul(h));
      const neg = c.mul(e).mul(g).add(a.mul(f).mul(h)).add(b.mul(d).mul(i));
      return {
        metodo: "Regla de Sarrus",
        valor: pos.sub(neg),
        expresion: `(${pos}) − (${neg})`,
      };
    }
    return null;
  }

  function operacionFila(i, factor, fila) {
    const neg = factor.signo() < 0;
    const k = factor.abs();
    const termino = k.esUno() ? `F${fila + 1}` : `(${k})·F${fila + 1}`;
    return `F${i + 1} ← F${i + 1} ${neg ? "+" : "−"} ${termino}`;
  }

  // ------------------------------------------------------------------
  // Sistemas de ecuaciones lineales (m ecuaciones, n variables, m <= n)
  // ------------------------------------------------------------------
  const INCOMPATIBLE = "incompatible";
  const DETERMINADO = "compatible determinado";
  const INDETERMINADO = "compatible indeterminado";

  function validarSistema(a) {
    const m = a.length;
    if (m === 0) throw new Error("La matriz no puede estar vacía.");
    const cols = a[0].length;
    if (cols < 2) throw new Error("La matriz debe tener al menos 2 columnas (una variable y el término independiente).");
    a.forEach((f, i) => {
      if (f.length !== cols) throw new Error(`La fila ${i + 1} tiene ${f.length} elementos; se esperaban ${cols}.`);
    });
    const n = cols - 1;
    if (m > n) throw new Error(`El número de ecuaciones (m = ${m}) debe ser menor o igual que el número de variables (n = ${n}).`);
    return [m, n];
  }

  /**
   * Resuelve el sistema dado por su matriz ampliada m x (n + 1).
   * Fase 1: eliminación hacia adelante (forma escalonada).
   * Fase 2: sustitución hacia atrás, con parámetros para las variables libres.
   */
  function resolverSistema(entrada, opciones = {}) {
    const pivoteoParcial = opciones.pivoteoParcial !== false;
    const original = aMatrizQ(entrada);
    const [m, n] = validarSistema(original);
    const a = copiar(original);
    const pasos = [];
    const pivotes = [];
    let fila = 0;

    for (let col = 0; col < n && fila < m; col++) {
      const p = elegirPivote(a, fila, m, col, pivoteoParcial);
      if (p < 0) continue; // columna sin pivote: variable libre
      if (p !== fila) {
        [a[fila], a[p]] = [a[p], a[fila]];
        pasos.push({ op: `F${fila + 1} ↔ F${p + 1}`, nota: "Intercambio de filas.", m: copiar(a), piv: [fila, col], filas: [fila, p], tipo: "swap" });
      }
      for (let i = fila + 1; i < m; i++) {
        if (a[i][col].esCero()) continue;
        const f = a[i][col].div(a[fila][col]);
        for (let j = col; j <= n; j++) a[i][j] = a[i][j].sub(f.mul(a[fila][j]));
        pasos.push({ op: operacionFila(i, f, fila), nota: `Se anula el elemento de la columna x${col + 1}.`, m: copiar(a), piv: [fila, col], filas: [i] });
      }
      pivotes.push(col);
      fila++;
    }

    const rango = pivotes.length;
    let rangoAmpliada = rango;
    for (let i = rango; i < m; i++) if (!a[i][n].esCero()) { rangoAmpliada = rango + 1; break; }

    const res = { m, n, rango, rangoAmpliada, escalonada: a, pivotes, pasos, original, solucion: null, libres: [] };
    if (rangoAmpliada !== rango) { res.tipo = INCOMPATIBLE; return res; }

    // Sustitución hacia atrás. Cada variable es una expresión afín:
    // { cte: Q, coef: Map(índice variable libre -> Q) }
    const libres = [];
    for (let j = 0; j < n; j++) if (!pivotes.includes(j)) libres.push(j);
    const sol = new Array(n);
    for (const j of libres) sol[j] = { cte: CERO, coef: new Map([[j, UNO]]) };

    for (let i = rango - 1; i >= 0; i--) {
      const c = pivotes[i];
      const piv = a[i][c];
      let cte = a[i][n].div(piv);
      const coef = new Map();
      for (let j = c + 1; j < n; j++) {
        if (a[i][j].esCero()) continue;
        const k = a[i][j].div(piv).neg();
        cte = cte.add(k.mul(sol[j].cte));
        for (const [t, v] of sol[j].coef) coef.set(t, (coef.get(t) || CERO).add(k.mul(v)));
      }
      for (const [t, v] of [...coef]) if (v.esCero()) coef.delete(t);
      sol[c] = { cte, coef };
    }

    res.solucion = sol;
    res.libres = libres;
    res.tipo = rango === n ? DETERMINADO : INDETERMINADO;
    return res;
  }

  /** Texto de una expresión afín, p. ej. "10 − 2·t1 − 3·t2". */
  function formatoExpresion(expr, nombres) {
    const partes = [];
    for (const j of [...expr.coef.keys()].sort((x, y) => x - y)) {
      const v = expr.coef.get(j);
      const k = v.abs();
      partes.push([v.signo() < 0 ? "−" : "+", k.esUno() ? nombres[j] : `${k}·${nombres[j]}`]);
    }
    let texto;
    if (!expr.cte.esCero() || partes.length === 0) {
      texto = expr.cte.signo() < 0 ? `−${expr.cte.abs()}` : `${expr.cte}`;
    } else {
      const [s, t] = partes.shift();
      texto = s === "+" ? t : `−${t}`;
    }
    for (const [s, t] of partes) texto += ` ${s} ${t}`;
    return texto;
  }

  /** Lista [{variable, valor, libre}] lista para mostrar. */
  function formatoSolucion(res) {
    if (!res.solucion) return null;
    const nombres = {};
    res.libres.forEach((j, k) => { nombres[j] = `t${k + 1}`; });
    return res.solucion.map((expr, j) => ({
      variable: `x${j + 1}`,
      valor: res.libres.includes(j) ? nombres[j] : formatoExpresion(expr, nombres),
      libre: res.libres.includes(j),
    }));
  }

  /** Evalúa la solución con valores dados para los parámetros (para comprobar). */
  function evaluarSolucion(res, valoresParametros) {
    const params = new Map(res.libres.map((j, k) => [j, parsear(valoresParametros[k] ?? 0)]));
    return res.solucion.map((e) => {
      let v = e.cte;
      for (const [j, c] of e.coef) v = v.add(c.mul(params.get(j)));
      return v;
    });
  }

  global.Gauss = {
    Q, CERO, UNO, parsear, aMatrizQ, determinante, comprobacionDirecta,
    resolverSistema, formatoSolucion, formatoExpresion, evaluarSolucion,
    INCOMPATIBLE, DETERMINADO, INDETERMINADO,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = global.Gauss;
})(typeof globalThis !== "undefined" ? globalThis : this);
