(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const out = [];
  for (let e = document.getElementById("ctx-below-1"); e; e = e.parentElement) {
    const c = getComputedStyle(e), r = e.getBoundingClientRect();
    out.push([e.tagName, String(e.className).slice(0, 40), Math.round(r.left), Math.round(r.right),
              e.scrollWidth, e.clientWidth, c.overflowX, c.width, e.style.width]);
  }
  const far = [];
  for (const e of document.querySelectorAll("*")) {
    const r = e.getBoundingClientRect();
    if (r.right > 1400 || r.left > 1400) far.push([e.tagName, String(e.className.baseVal === undefined ? e.className : e.className.baseVal).slice(0, 30), Math.round(r.left), Math.round(r.right), Math.round(r.width)]);
  }
  return {chain: out, far: far.slice(0, 20), scroll: document.documentElement.scrollWidth,
          body: document.body.scrollWidth};
})()
