from loguru import logger
from fastapi import WebSocket
import asyncio


class CallbackHandler:
    def __init__(self, websocket: WebSocket):
        self.websocket = websocket

    async def send(self, assistant_message):
        send = self.websocket.send_json
        if assistant_message.type == "UserMessage":
            message = f"User: {assistant_message.content}"
            logger.info(message)
            await send({"type": "response", "message": message})
        elif assistant_message.type == "AssistantMessage":

            if assistant_message.thought is not None:
                _str = f"Model thought: {assistant_message.thought}"
                logger.info(_str)
                await send(
                    {
                        "type": "response",
                        "message": _str,
                    }
                )
            if isinstance(assistant_message.content, list):
                for item in assistant_message.content:
                    if hasattr(item, "name") and hasattr(item, "arguments"):
                        _str = f"Function call: {item.name} with args {item.arguments}"
                        logger.info(_str)
                        msg = {"type": "response", "message": _str}
                        if "smiles" in item.arguments:
                            msg["smiles"] = item.arguments["smiles"]
                        await send(msg)

                    else:

                        logger.info(f"Model: {item}")
        elif assistant_message.type == "FunctionExecutionResultMessage":

            for result in assistant_message.content:
                if result.is_error:
                    message = (
                        f"Function {result.name} errored with output: {result.content}"
                    )
                    logger.error(message)
                else:
                    message = f"Function {result.name} returned: {result.content}"
                    logger.info(message)
        else:
            message = f"Model: {assistant_message.message.content}"
            logger.info(message)

    def __call__(self, assistant_message):
        asyncio.create_task(self.send(assistant_message))
