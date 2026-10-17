(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const fig = document.querySelector('#figure-8'), d = fig.querySelector('.js-plotly-plot');
  const form = fig.querySelector('form.view');
  const read = () => Array.from(d.querySelectorAll('.legend .traces')).slice(0, 4).map(g => [
    (g.querySelector('.legendtext') || {}).textContent,
    Array.from(g.querySelectorAll('path')).map(p => p.getAttribute('style')).join(' | ')]);
  const out = {build: read()};
  for (const v of ['1', '0', '2']) {
    form.elements.view.value = v; form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
    await wait(400);
    out['view' + v] = read();
  }
  return out;
})()
