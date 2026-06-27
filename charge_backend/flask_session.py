###############################################################################
## Copyright 2025-2026 Lawrence Livermore National Security, LLC.
## See the top-level LICENSE file for details.
##
## SPDX-License-Identifier: Apache-2.0
###############################################################################
"""
Definition of a FLASK user session object.
"""

import argparse
from typing import Any
import uuid
from fastapi import WebSocket
from loguru import logger

from lc_conductor.local_mcp_proxy import resolve_local_mcp_response
from lc_conductor.session import UserSession

from charge_backend.backend_manager import FlaskActionManager
from charge_backend import prompt_debugger


class FlaskUserSession(UserSession):
    """
    A representation of a persistent FLASK Copilot user session. Includes the
    action manager, which in turn includes information about the currently-running
    experiment, UI state, and AI agents and their context.
    """

    action_manager: FlaskActionManager

    def __init__(
        self, username: str, args: argparse.Namespace, websocket: WebSocket
    ) -> None:
        super().__init__(username, str(uuid.uuid4()), websocket, action_manager=None)

        action_manager = FlaskActionManager(
            self.websocket,
            args,
            username,
        )
        self.action_manager: FlaskActionManager = action_manager

    async def event_loop(self):
        # Before the event loop starts, synchronize the orchestrator model config
        await self.action_manager.report_orchestrator_config()

        # Run event loop normally
        return await super().event_loop()

    async def handle_action(self, action: str, data: dict[str, Any]):
        """
        Handles a received message from the FLASK Copilot Web UI.
        First we look for prompt breakpoint and local MCP proxy responses,
        because they are received on an asynchronous thread during an existing
        request. Then, we try to dispatch from the FLASK Action Manager.

        :param action: The name of the action to perform.
        :param data: The JSON data received with the message.
        """
        if action == "prompt-breakpoint-response":  # AI debugging
            prompt_debugger.DEBUG_PROMPT_RESPONSES[self.websocket].set_result(data)
            return

        if action == "local-mcp-response":
            if not resolve_local_mcp_response(self.websocket, data):
                logger.warning(f"Received unmatched local MCP response: {data}")
            return

        if not self.action_manager.has_handler(action):
            logger.warning(f"Unknown action received: {action} with data {data}")

        return await self.action_manager.dispatch(action, data)
