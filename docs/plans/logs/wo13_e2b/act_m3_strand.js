(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const form = document.querySelector('form.strand'), fig = form.closest('.fig');
  const d = fig.querySelector('.js-plotly-plot');
  const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
  const read = () => {
    const out = {summary: fig.querySelector('.strand-summary').textContent, calls: {},
                 moved: !!fig.querySelector('a[data-export][role=button]'),
                 lines: cfg.lines.map(i => d.data[i].y[0])};
    cfg.calls.forEach(c => {
      out.calls[c] = (d.data[cfg.traces[c]].customdata || []).map(r => r[0]);
    });
    return out;
  };
  const out = {loaded: read(), settings: []};
  for (const [st, un] of [['0.6', '0.3'], ['0.95', '0.05'], ['0.8', '0.1']]) {
    form.elements.stranded.value = st; form.elements.unstranded.value = un;
    form.elements.unstranded.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(150); out.settings.push(read());
  }
  out.box = document.documentElement.clientWidth; out.scroll = document.documentElement.scrollWidth;
  return out;
})()