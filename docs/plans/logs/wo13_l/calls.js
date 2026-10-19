(async () => {
  /* Every control bar (.calls, .survival, .brush) at this width: a label, select text or button
     cut by its own box or reaching past the bar. A select's shown text is measured with the
     select's own font against its content box less the arrow. */
  await new Promise(r => setTimeout(r, 4000));
  const ctx = document.createElement("canvas").getContext("2d"), out = [];
  for (const bar of document.querySelectorAll(".calls, .survival, .brush")) {
    if (!bar.getClientRects().length) continue;
    const box = bar.getBoundingClientRect();
    for (const e of bar.querySelectorAll("label, select, input, button, legend, span")) {
      if (!e.getClientRects().length) continue;
      const r = e.getBoundingClientRect(), c = getComputedStyle(e);
      let cut = r.right > box.right + 0.5 || r.left < box.left - 0.5;
      let shown = (e.textContent || "").trim();
      if (e.tagName === "SELECT" && e.selectedIndex >= 0) {
        shown = e.options[e.selectedIndex].text;
        ctx.font = `${c.fontWeight} ${c.fontSize} ${c.fontFamily}`;
        const room = e.clientWidth - parseFloat(c.paddingLeft) - parseFloat(c.paddingRight) - 16;
        cut = cut || ctx.measureText(shown).width > room;
      } else if (e.tagName !== "INPUT") {
        cut = cut || e.scrollWidth > e.clientWidth + 1;
      }
      if (cut) out.push([String(bar.className).slice(0, 30), e.tagName, Math.round(r.left),
                         Math.round(r.right), Math.round(box.right), shown.slice(0, 60)]);
    }
  }
  return out;
})()
