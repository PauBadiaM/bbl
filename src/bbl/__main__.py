"""``python -m bbl`` -- the same entry point as the ``bbl`` console script."""

from .cli import main

raise SystemExit(main())
