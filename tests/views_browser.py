"""Walk every option of every ``view`` select on a page in headless Chromium (D-073, WO6 to WO10).

A plain view sets its select to each option; a grid (``ctx.report.lookups``) sets its selects to
each combination it names, through change events as a reader's would fire. For each option it
reads what plotly drew (trace names, x, y, annotation texts), the title, lead, caption and live
line, and whether the PDF and SVG buttons say the view is moved; then it sets the build back.
It also measures the page box against the window, the theme drawn and console errors.
"""

from pathlib import Path

from test_export_browser import _browse

WALK = """(async () => {
  const wait = ms => new Promise(r => setTimeout(r, ms));
  await wait(1500);
  const forms = Array.from(document.querySelectorAll('form.view'));
  const out = [];
  for (const form of forms) {
    const fig = form.closest('.fig'), d = fig.querySelector('.js-plotly-plot');
    const cfg = JSON.parse(document.getElementById(form.getAttribute('data-controls')).textContent);
    const words = s => (fig.querySelector(s) || {}).textContent || '';
    const read = () => ({
      traces: d.data.map(t => ({name: t.name, x: t.x, y: t.y, colour: t.marker && t.marker.color})),
      notes: ((d.layout || {}).annotations || []).map(a => a.text),
      hidden: ((d.layout || {}).annotations || []).filter(a => a.visible === false).length,
      ticks: ((d.layout || {}).xaxis || {}).ticktext || null,
      title: words('.fig-title'), lead: words('.fig-lead'), caption: words('figcaption'),
      line: words('.view-summary'), moved: !!fig.querySelector('a[data-export][role=button]'),
      wide: d.getBoundingClientRect().width});
    const selects = Array.from(form.querySelectorAll('select[data-axis]'));
    const loaded = read();
    const set = async (i) => {
      if (selects.length) {
        const key = Object.keys(cfg.grid).find(k => cfg.grid[k] === i).split(',');
        selects.forEach((s, j) => { s.value = key[j]; });
        selects[0].dispatchEvent(new Event('change', {bubbles: true}));
      } else {
        form.elements.view.value = String(i);
        form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await wait(60);
    };
    const options = [];
    for (let i = 0; i < cfg.options.length; i++) { await set(i); options.push(read()); }
    await set(0);
    out.push({figure: fig.id, label: cfg.label, grid: !!selects.length, options: options,
              loaded: loaded, back: read(), built: cfg.options.map(o => o.figure.data.map(t => ({
                name: t.name, x: t.x, y: t.y})))});
  }
  return {views: out, theme: document.documentElement.getAttribute('data-theme'),
          box: document.documentElement.clientWidth,
          scroll: document.documentElement.scrollWidth};
})()"""

SIZES = (("wide_light", 1280, "light"), ("phone_dark", 375, "dark"))


def walk(page: Path, downloads: Path, sizes=SIZES) -> dict[str, dict]:
    """Every view's options on ``page`` at each (name, width, scheme); results by name."""
    tasks = [{"name": name, "url": page.as_uri(), "width": width, "height": 900,
              "scheme": scheme, "act": WALK} for name, width, scheme in sizes]  # fmt: skip
    return _browse(tasks, downloads)


def exported(page: Path, downloads: Path, figure: str, option: int) -> dict:
    """Set one view on ``page`` to ``option`` and press its PDF and SVG buttons."""
    act = """(async () => {
      const wait = ms => new Promise(r => setTimeout(r, ms));
      await wait(1500);
      const form = document.querySelector('FIGURE form.view'), at = OPTION;
      const id = form.getAttribute('data-controls');
      const cfg = JSON.parse(document.getElementById(id).textContent);
      const selects = Array.from(form.querySelectorAll('select[data-axis]'));
      if (selects.length) {
        const key = Object.keys(cfg.grid).find(k => cfg.grid[k] === at).split(',');
        selects.forEach((s, j) => { s.value = key[j]; });
        selects[0].dispatchEvent(new Event('change', {bubbles: true}));
      } else {
        form.elements.view.value = String(at);
        form.elements.view.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await wait(200);
      return cfg.options[at].export ? cfg.options[at].export.caption : null;
    })()""".replace("FIGURE", figure).replace("OPTION", str(int(option)))
    (result,) = _browse([{"name": "export", "url": page.as_uri(), "width": 1280, "act": act,
                          "figure": figure}], downloads).values()  # fmt: skip
    return result
