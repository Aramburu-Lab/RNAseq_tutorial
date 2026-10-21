const fs = require("fs");
eval(fs.readFileSync(process.argv[2], "utf8") + ";globalThis.ctxUnpack=ctxUnpack;");
const p = fs.readFileSync(process.argv[3], "utf8");
for (let k = 0; k < 5; k++) {
  const t = performance.now();
  const text = atob(p), bytes = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i++) bytes[i] = text.charCodeAt(i);
  const t1 = performance.now();
  const s = new TextDecoder().decode(ctxUnpack.inflate(bytes));
  console.log((t1 - t).toFixed(1), "ms b64;", (performance.now() - t1).toFixed(1), "ms inflate+decode", s.length);
}
