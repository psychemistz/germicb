"""
Markdown report generation helpers.
"""

import time
from typing import List, Optional


def markdown_table(headers: List[str], rows: List[List[str]]) -> str:
    """Generate a markdown table string."""
    lines = []
    lines.append('| ' + ' | '.join(headers) + ' |')
    lines.append('|' + '|'.join('------' for _ in headers) + '|')
    for row in rows:
        lines.append('| ' + ' | '.join(str(c) for c in row) + ' |')
    return '\n'.join(lines)


def report_header(title: str, subtitle: Optional[str] = None) -> str:
    """Generate a standard report header with title and timestamp."""
    header = f"# {title}\n"
    if subtitle:
        header += f"## {subtitle}\n"
    header += f"\n_Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}_\n\n"
    return header
