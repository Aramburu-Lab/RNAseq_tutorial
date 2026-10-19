(async () => {
  await new Promise(r => setTimeout(r, 4000));
  const w = document.documentElement.clientWidth, far = [];
  const walk = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
  for (let n = walk.nextNode(); n; n = walk.nextNode()) {
    if (!n.textContent.trim()) continue;
    const range = document.createRange(); range.selectNodeContents(n);
    const r = range.getBoundingClientRect();
    if (r.right > w + 0.5) {
      const p = n.parentElement;
      far.push([p.tagName, String(p.className.baseVal === undefined ? p.className : p.className.baseVal).slice(0, 30),
                Math.round(r.left), Math.round(r.right), n.textContent.trim().slice(0, 60)]);
    }
  }
  /* Elements whose own scroll width passes their box, outermost first, outside sideways scrollers. */
  const spill = [];
  for (const e of document.querySelectorAll("body *")) {
    const c = getComputedStyle(e);
    if (c.overflowX !== "visible") continue;
    if (e.scrollWidth > e.clientWidth + 1 && e.clientWidth > 0) {
      spill.push([e.tagName, String(e.className.baseVal === undefined ? e.className : e.className.baseVal).slice(0, 30),
                  e.clientWidth, e.scrollWidth, (e.textContent || "").trim().slice(0, 40)]);
    }
  }
  return {far: far.slice(0, 15), spill: spill.slice(-15)};
})()
