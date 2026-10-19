(async () => {
  /* Every plot's tick labels and axis titles cut by the plot's own svg (the left or right
     edge of the drawn figure), each figure shown first by its tab where it sits in one. */
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(3000);
  const out = [];
  for (const fig of document.querySelectorAll("figure.fig")) {
    if (fig.hidden) {
      const tab = document.querySelector(`[role=tab][aria-controls="${fig.id}"]`);
      if (tab) { tab.click(); await wait(1200); }
    }
    for (const d of fig.querySelectorAll(".js-plotly-plot")) {
      const svg = d.querySelector("svg.main-svg");
      if (!svg) continue;
      const box = svg.getBoundingClientRect();
      for (const t of d.querySelectorAll(".ytick text, .xtick text, .g-ytitle text, .g-xtitle text, .legendtext")) {
        const r = t.getBoundingClientRect();
        if (r.width && (r.left < box.left - 0.5 || r.right > box.right + 0.5)) {
          out.push([fig.querySelector(".fig-title").textContent.slice(0, 50), t.textContent.slice(0, 50),
                    Math.round(r.left - box.left), Math.round(r.right - box.left), Math.round(box.width)]);
        }
      }
    }
  }
  return out;
})()
