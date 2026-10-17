(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const fig = document.querySelector('#figure-8'), d = fig.querySelector('.js-plotly-plot');
  const form = fig.querySelector('form.view');
  const read = () => d._fullData.slice(0, 4).map(t => [t.name, t.line && t.line.color, t.line && t.line.dash,
    t.contours && t.contours.start, t._input && t._input.line && t._input.line.color]);
  const out = {build: read()};
  form.elements.view.value = '1'; form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
  await wait(400);
  out.other = read();
  out.paths = Array.from(d.querySelectorAll('.contour path.contourlevel, .contour path')).slice(0, 6)
    .map(p => p.getAttribute('style'));
  return out;
})()
