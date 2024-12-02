from typing import Optional
from pathlib import Path


def parse_filter(filter: Path | list[str] | set[str] | None = None) -> Optional[set[str]]:
    """Converts the input into a set of strings to use as a filter.

    Args:
        filter (Path | list[str] | set[str] | None, optional):
            The filter can be:
            - A `Path` to a file where each line contains a string.
            - A `list` or `set` of strings.
            - `None`, in which case the function will return `None`.

    Returns:
        Optional[set[str]]:
            A set of strings to use as a filter, or `None` if the input is `None`.
    """
    if not filter:
        return None

    if isinstance(filter, Path):
        return set(filter.read_text().strip().split("\n"))

    return set(filter)


def format_fig(fig):
    """Formats a plotly figure in a nicer way."""
    fig.update_layout(
        width=1200,
        height=550,
        plot_bgcolor="white",
        title_font_color="black",
        font=dict(
            family="Linux Libertine Display O",
            size=18,
            color="black",
        ),
    )
    gridcolor = "#dddddd"
    fig.update_xaxes(gridcolor=gridcolor)
    fig.update_yaxes(gridcolor=gridcolor)

    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True, ticks='outside')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True, ticks='outside')
