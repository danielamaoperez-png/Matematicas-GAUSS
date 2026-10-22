// Service worker: guarda la app en el teléfono para que funcione sin conexión.
// Cambia VERSION cada vez que se publique una versión nueva de la app.
const VERSION = "gauss-v1";
const ARCHIVOS = [
  "./",
  "./index.html",
  "./gauss.js",
  "./manifest.webmanifest",
  "./icons/icon-180.png",
  "./icons/icon-192.png",
  "./icons/icon-512.png",
];

self.addEventListener("install", (e) => {
  e.waitUntil(caches.open(VERSION).then((c) => c.addAll(ARCHIVOS)).then(() => self.skipWaiting()));
});

self.addEventListener("activate", (e) => {
  e.waitUntil(
    caches.keys()
      .then((claves) => Promise.all(claves.filter((k) => k !== VERSION).map((k) => caches.delete(k))))
      .then(() => self.clients.claim())
  );
});

// Primero responde desde la caché (rápido y sin conexión) y actualiza en segundo plano.
self.addEventListener("fetch", (e) => {
  if (e.request.method !== "GET") return;
  e.respondWith(
    caches.open(VERSION).then((cache) =>
      cache.match(e.request, { ignoreSearch: true }).then((guardada) => {
        const red = fetch(e.request)
          .then((resp) => {
            if (resp.ok || resp.type === "opaque") cache.put(e.request, resp.clone());
            return resp;
          })
          .catch(() => guardada);
        return guardada || red;
      })
    )
  );
});
