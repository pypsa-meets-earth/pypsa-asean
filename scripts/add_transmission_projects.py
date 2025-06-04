# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Add transmission projects to the network.
"""

from pathlib import Path

import pandas as pd
import pypsa

from _helpers import configure_logging, create_logger  # make sure this is still present

logger = create_logger(__name__)


def attach_transmission_projects(
    n: pypsa.Network, transmission_projects: list[str]
) -> None:
    logger.info("Adding transmission projects to network.")
    for path in transmission_projects:
        path = Path(path)
        df = pd.read_csv(path, index_col=0, dtype={"bus0": str, "bus1": str})
        if df.empty:
            continue
        if "new_buses" in path.name:
            n.add("Bus", df.index, **df)
        elif "new_lines" in path.name:
            n.add("Line", df.index, **df)
        elif "new_links" in path.name:
            n.add("Link", df.index, **df)
        elif "adjust_lines" in path.name:
            n.lines.update(df)
        elif "adjust_links" in path.name:
            n.links.update(df)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("add_transmission_projects_and_dlr")
    configure_logging(snakemake)  # pylint: disable=E0606

    n = pypsa.Network(snakemake.input.network)

    params = snakemake.params

    if params["transmission_projects"]["enable"]:
        attach_transmission_projects(n, snakemake.input.transmission_projects)

    n.export_to_netcdf(snakemake.output[0])
