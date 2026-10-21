# ruff: noqa: D103, E501, PLR0911, PLR2004, SIM115, B007 - a proof record, run once and kept as run
"""wo13-z1: what every byte of an M10 page holds, read-only.

Run: python measure.py <page.html>. Every byte of the page lands in exactly one row: each script
and style by what it carries, each figure split into its div and script wrapper, its traces,
its layout (the plotly template apart) and its config, its rows and its controls (the survival
panel's by key), and the HTML outside every script and style as text. The rows add up to the
page's size, and the script says so. Beside each part: how it is encoded, and its size packed
as builder.pack packs (gzip with no timestamp, then base64), so a raw part's saving is read off.
Duplication: the plotly template every figure carries beside ctx-templates, the patient barcodes
and gene symbols counted in each part that holds them.
"""

import base64
import gzip
import json
import re
import sys
from collections import Counter

PART = re.compile(r"<(script|style)([^>]*)>(.*?)</\1>", re.S)
NEWPLOT = re.compile(r"Plotly\.newPlot\(\s*\"([^\"]+)\",\s*")
BARCODE = re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")


def packed(text: str) -> int:
    """Bytes of ``text`` gzipped with no timestamp and in base64, as builder.pack writes it."""
    return len(base64.b64encode(gzip.compress(text.encode(), mtime=0)))


def size(text: str) -> int:
    return len(text.encode())


def figure_rows(name: str, body: str) -> list[tuple[str, str, int, int]]:
    """A figure script as its wrapper, traces, layout without its template, template, config."""
    match = NEWPLOT.search(body)
    if match is None:
        return [(name, "inline JS", size(body), packed(body))]
    held = re.compile(r'\.\.\.ctxUnpack\("([^"]*)"\)').match(body, match.end())
    if held:  # packed by figures.pack_json (Z.1)
        return [
            (f"{name} script wrapper", "inline JS", size(body) - len(held.group(1)), 0),
            (f"{name} three arguments", "DEFLATE+base64 JSON", len(held.group(1)), 0),
        ]
    decoder, at, values, spans = json.JSONDecoder(), match.end(), [], []
    for _ in range(3):
        start = at
        value, at = decoder.raw_decode(body, at)
        values.append(value)
        spans.append(body[start:at])
        at = re.compile(r"\s*,?\s*").match(body, at).end()
    data, layout, _ = values
    template = json.dumps(layout.pop("template", {}), separators=(",", ":"))
    wrapper = size(body) - sum(size(s) for s in spans)
    rest = size(spans[1]) - size(template)
    args = ",".join(spans)
    return [
        (f"{name} script wrapper", "inline JS", wrapper, 0),
        (f"{name} traces ({len(data)})", "raw JSON", size(spans[0]), 0),
        (f"{name} layout, template apart", "raw JSON", rest, 0),
        (f"{name} layout template", "raw JSON", size(template), 0),
        (f"{name} config; packed: all three arguments", "raw JSON", size(spans[2]), packed(args)),
    ]


def controls_rows(name: str, body: str) -> list[tuple[str, str, int, int]]:
    """A survival panel's settings, key by key; any other controls as one row."""
    settings = json.loads(body)
    if settings.get("kind") != "survival":
        return [(name, "raw JSON", size(body), packed(body))]
    out, used = [], 0
    for key, value in settings.items():
        text = json.dumps(value)
        enc = "gzip+base64 int16 byte planes" if key == "block" else "raw JSON"
        out.append((f"{name} .{key}", enc, size(text), 0 if key == "block" else packed(text)))
        used += size(text)
    out.append((f"{name} JSON keys and separators", "raw JSON", size(body) - used, 0))
    return out


def classify(attrs: str, body: str, index: int) -> list[tuple[str, str, int, int]]:
    """The rows one script or style element's content makes."""
    ident = re.search(r'id="([^"]+)"', attrs)
    name = ident.group(1) if ident else f"script at {index}"
    if attrs == "style":
        return [("CSS", "raw CSS", size(body), packed(body))]
    if body.startswith("/**\n* plotly.js"):
        return [("plotly.js", "raw JS", size(body), packed(body))]
    if body.startswith("/*! jsPDF"):
        return [("export libraries (jsPDF, svg2pdf)", "raw JS", size(body), packed(body))]
    if "PLOTLYENV" in body[:200]:
        match = NEWPLOT.search(body)
        return figure_rows(match.group(1) if match else name, body)
    if name.startswith("ctx-controls-"):
        return controls_rows(name, body)
    if name.startswith("ctx-data-") and "octet-stream" in attrs:
        return [(f"{name} rows", "DEFLATE+base64 JSON", size(body), 0)]
    if name.startswith("ctx-data-"):
        return [(f"{name} rows", "raw JSON (columnar)", size(body), packed(body))]
    if name == "ctx-export":
        return [("ctx-export: fonts and settings", "base64 fonts in JSON", size(body), 0)]
    if ident:
        return [(name, "raw JSON", size(body), packed(body))]
    kind = "page script" if size(body) > 50_000 else "small inline script"
    return [(f"{kind} ({name})", "raw JS", size(body), packed(body))]


def main(path: str) -> None:
    html = open(path, encoding="utf-8").read()
    rows: list[tuple[str, str, int, int]] = []
    outside, last = 0, 0
    for m in PART.finditer(html):
        tag, attrs, body = m.group(1), m.group(2), m.group(3)
        outside += size(html[last : m.start(3)]) + size(f"</{tag}>")
        last = m.end()
        rows += classify("style" if tag == "style" else attrs, body, m.start())
    outside += size(html[last:])
    rows.append(
        ("HTML outside scripts and styles: text, tags, tags of the parts", "HTML", outside, 0)
    )
    total = sum(r[2] for r in rows)
    print(f"page {path}\nbytes {size(html):,}; rows sum {total:,}; equal: {total == size(html)}")
    print(f"{'part':<62} {'encoding':<32} {'bytes':>11} {'packed':>9}")
    for name, enc, n, p in rows:
        print(f"{name:<62} {enc:<32} {n:>11,} {p if p else '':>9}")
    kinds = Counter()
    for name, enc, n, _ in rows:
        kinds[enc] += n
    print("\nby encoding")
    for enc, n in kinds.most_common():
        print(f"  {enc:<32} {n:>11,}")
    duplication(html)


def duplication(html: str) -> None:
    """The template in each figure against ctx-templates; barcodes and symbols per part."""
    templates = json.loads(re.search(r'id="ctx-templates">(.*?)</script>', html, re.S).group(1))
    print("\nduplication")
    copies = 0
    for m in PART.finditer(html):
        if m.group(1) == "script" and "PLOTLYENV" in m.group(3)[:200]:
            body = m.group(3)
            nm = NEWPLOT.search(body)
            if body[nm.end() :].startswith("...ctxUnpack"):
                continue  # packed: its template is inside the packed text
            decoder = json.JSONDecoder()
            _, at = decoder.raw_decode(body, nm.end())
            at = re.compile(r"\s*,\s*").match(body, at).end()
            layout, _ = decoder.raw_decode(body, at)
            copies += size(json.dumps(layout.get("template", {}), separators=(",", ":")))
    themes = {k: size(json.dumps(v)) for k, v in templates.items()}
    print(
        f"  plotly template inside the figures, all copies: {copies:,} B; ctx-templates: {themes}"
    )
    for label, pattern in (("patient barcodes", BARCODE),):
        for m in PART.finditer(html):
            found = pattern.findall(m.group(3))
            if found:
                ident = re.search(r'id="([^"]+)"', m.group(2))
                nm = NEWPLOT.search(m.group(3)[:400])
                where = ident.group(1) if ident else (nm.group(1) if nm else "script")
                print(
                    f"  {label} in {where}: {len(found):,} occurrences, {len(set(found)):,} distinct"
                )
    survival = re.search(r'id="ctx-controls-\d+">(\{"kind": "survival".*?)</script>', html, re.S)
    if survival:
        cfg = json.loads(survival.group(1))
        listed = {row[2][0] for row in cfg["screen"]}
        print(
            f"  survival panel: {cfg['genes']:,} genes by {cfg['patients']:,} patients in the block;"
        )
        print(f"    symbols {len(cfg['symbols']):,}, screen rows {len(cfg['screen']):,}, ")
        print(f"    screen genes that are panel genes: {len(listed & set(cfg['symbols'])):,}")
        print(
            f"    symbols as their own list: {size(json.dumps(cfg['symbols'])):,} B, repeated in the screen rows"
        )


if __name__ == "__main__":
    main(sys.argv[1])
