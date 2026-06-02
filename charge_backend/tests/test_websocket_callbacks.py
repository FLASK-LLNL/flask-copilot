import asyncio

from charge_backend.backend_helper_funcs import CallbackHandler
from lc_conductor.callback_logger import CallbackLogger


class FakeWebSocket:
    def __init__(self):
        self.messages: list[dict] = []

    async def send_json(self, payload: dict) -> None:
        self.messages.append(payload)


class FakeAssistantMessage:
    type = "UserMessage"
    content = "hello"


def test_callback_logger_sends_immediately_and_in_order():
    async def run() -> None:
        websocket = FakeWebSocket()
        clogger = CallbackLogger(websocket, source="test")

        await clogger.info("first")
        await websocket.send_json({"type": "complete"})

        assert websocket.messages == [
            {
                "type": "response",
                "message": {"source": "test", "message": "first"},
            },
            {"type": "complete"},
        ]

    asyncio.run(run())


def test_callback_handler_drain_flushes_background_send():
    async def run() -> None:
        websocket = FakeWebSocket()
        callback = CallbackHandler(websocket, agent_key="test-agent")

        callback(FakeAssistantMessage())
        assert websocket.messages == []

        await callback.drain()

        assert websocket.messages == [
            {
                "type": "response",
                "message": {"message": "[test-agent] User: hello"},
            }
        ]

    asyncio.run(run())


def test_callback_handler_notifies_agent_update_without_tagging_sidebar_messages():
    async def run() -> None:
        websocket = FakeWebSocket()
        update_count = 0

        async def on_update() -> None:
            nonlocal update_count
            update_count += 1

        callback = CallbackHandler(
            websocket,
            agent_key="reaction:node_1",
            on_agent_update=on_update,
        )

        await callback.on_tool_result("lookup_reaction", {"ok": True})

        assert "agentKey" not in websocket.messages[0]
        assert websocket.messages[0]["message"]["source"] == "lookup_reaction"
        assert update_count == 1

    asyncio.run(run())
