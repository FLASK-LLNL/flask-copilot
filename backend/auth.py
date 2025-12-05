################################################################################
## Copyright 2025 Lawrence Livermore National Security, LLC. and Binghamton University.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
################################################################################

from fastapi import Header, HTTPException, status
from typing import Optional


async def get_forwarded_user(x_forwarded_user: Optional[str] = Header(None, alias="X-Forwarded-User")) -> str:
    """Extract authenticated user from X-Forwarded-User header"""
    if not x_forwarded_user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing X-Forwarded-User header. Authentication required.",
        )
    return x_forwarded_user
