"""Resolve the external S-LDSC resource/result root."""

from __future__ import annotations

import os
from pathlib import Path


def work_root() -> Path:
    configured = os.environ.get("SC_PCQTL_SLDSC_ROOT", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(__file__).resolve().parents[1]
