// Pruebas de la lógica de la app (app/gauss.js). Ejecutar: node --test tests/
const test = require("node:test");
const assert = require("node:assert/strict");
const G = require("../app/gauss.js");

const s = (q) => q.toString();

test("parsea enteros, decimales, comas, fracciones y signo menos tipográfico", () => {
  assert.equal(s(G.parsear("3")), "3");
  assert.equal(s(G.parsear("-2.5")), "-5/2");
  assert.equal(s(G.parsear("2,5")), "5/2");
  assert.equal(s(G.parsear("1/3")), "1/3");
  assert.equal(s(G.parsear("−4")), "-4");
  assert.equal(s(G.parsear("-.5")), "-1/2");
  assert.equal(s(G.parsear("")), "0");
  assert.equal(G.parsear("abc"), null);
  assert.equal(G.parsear("1/0"), null);
});

test("determinante 3x3 coincide con Sarrus", () => {
  const A = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]];
  const r = G.determinante(A);
  assert.equal(s(r.det), "-1");
  assert.equal(s(G.comprobacionDirecta(r.original).valor), "-1");
});

test("determinante con intercambio obligatorio y sin pivoteo parcial", () => {
  const r = G.determinante([[0, 1], [1, 0]], { pivoteoParcial: false });
  assert.equal(s(r.det), "-1");
  assert.equal(r.intercambios, 1);
});

test("determinante 4x4 y con fracciones", () => {
  const A = [[1, 2, 3, 4], [5, 6, 7, 8], [2, 6, 4, 8], [3, 1, 1, 2]];
  assert.equal(s(G.determinante(A).det), "72");
  assert.equal(s(G.determinante(A, { pivoteoParcial: false }).det), "72");
  assert.equal(s(G.determinante([["1/2", "1/3"], ["1/4", "1/5"]]).det), "1/60");
});

test("matriz singular", () => {
  const r = G.determinante([[1, 2, 3], [2, 4, 6], [1, 0, 1]]);
  assert.equal(s(r.det), "0");
  assert.ok(r.singular);
});

test("rechaza matrices no cuadradas y valores no numéricos", () => {
  assert.throws(() => G.determinante([[1, 2, 3], [4, 5, 6]]));
  assert.throws(() => G.determinante([["x", 1], [1, 1]]), /no es un número/);
});

function comprobar(A, res, params) {
  const x = G.evaluarSolucion(res, params);
  const M = G.aMatrizQ(A);
  for (const f of M) {
    let izq = G.CERO;
    for (let j = 0; j < x.length; j++) izq = izq.add(f[j].mul(x[j]));
    assert.ok(izq.eq(f[f.length - 1]));
  }
}

test("sistema compatible determinado", () => {
  const A = [[2, 1, -1, 8], [-3, -1, 2, -11], [-2, 1, 2, -3]];
  const r = G.resolverSistema(A);
  assert.equal(r.tipo, G.DETERMINADO);
  assert.deepEqual(G.formatoSolucion(r).map((v) => v.valor), ["2", "3", "−1"]);
});

test("sistema indeterminado con menos ecuaciones que variables", () => {
  const A = [[1, 2, -1, 1, 3], [2, 4, 1, -1, 3]];
  const r = G.resolverSistema(A);
  assert.equal(r.tipo, G.INDETERMINADO);
  assert.equal(r.libres.length, 2);
  for (const p of [[0, 0], [1, -2], ["7/3", 5]]) comprobar(A, r, p);
  const una = G.resolverSistema([[1, 2, 3, 4, 10]]);
  assert.equal(G.formatoSolucion(una)[0].valor, "10 − 2·t1 − 3·t2 − 4·t3");
});

test("sistema incompatible y validación m <= n", () => {
  const r = G.resolverSistema([[1, 1, 2], [1, 1, 3]]);
  assert.equal(r.tipo, G.INCOMPATIBLE);
  assert.equal(r.rango, 1);
  assert.equal(r.rangoAmpliada, 2);
  assert.throws(() => G.resolverSistema([[1, 2], [3, 4]]), /menor o igual/);
});
