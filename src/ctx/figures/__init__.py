"""Reusable figures (Phase 2D, D-020): a stable ID, a parameter model and ``ctx figure``.

``kinds`` says which figures have a spec and what each reads; ``record`` builds a figure's
``FigureSpec`` when a module writes it; ``panel`` turns that spec into the Reproduce panel;
``command`` is ``ctx figure``, which redraws a figure from its inputs alone, in the reader's
own folder. The ``draw_*`` modules call the module's own figure function on tables read back,
so a figure has one definition however it is drawn.
"""
