################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from datetime import datetime, timedelta


def next_last_modified(current: datetime | None) -> datetime:
    """Return a monotonic UTC timestamp compatible with second-precision stores."""
    now = datetime.utcnow().replace(microsecond=0)
    if current is None:
        return now

    current_value = current.replace(tzinfo=None, microsecond=0)
    if now <= current_value:
        return current_value + timedelta(seconds=1)
    return now
