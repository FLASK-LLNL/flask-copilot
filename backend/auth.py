################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from fastapi import Header
from typing import Optional
import os


async def get_forwarded_user(x_forwarded_user: Optional[str] = Header(None, alias="X-Forwarded-User")) -> str:
    """
    Extract authenticated user from X-Forwarded-User header.
    Falls back to 'default_user' for local development.
    """
    if not x_forwarded_user:
        # For local development, use a default user
        return os.getenv("DEFAULT_USER", "default_user")
    return x_forwarded_user
