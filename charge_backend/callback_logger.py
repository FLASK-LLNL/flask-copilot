import sys
from fastapi import WebSocket
from loguru import logger
from typing import Optional


# Define the callback function - will send message to the websocket if it is provided
async def handle_callback_log(message):
    record = message.record
    websocket = record["extra"].get("websocket", None)
    smiles = record["extra"].get("smiles", None)
    source = record["extra"].get("source", None)
    kwargs = {}
    if smiles:
        kwargs["smiles"] = smiles
    if websocket:
        # Timestamp is already included in the GUI window
        # timestamp = record["time"].isoformat(" ", timespec='seconds')
        msg = record["message"]
        level = record["level"].name
        if not source:
            LEVELS = {
               'DEBUG': 'Debug',
               'VERBOSE': 'Verbose',
               'INFO': 'Info',
               'WARN': 'Warning',
               'WARNING': 'Warning',
                'ERROR': 'Error',
               'EXCEPTION': 'Exception'
            }
            level_str = LEVELS.get(level, level)
            source = f"Logger ({level_str})"
        await websocket.send_json(
            {
                "type": "response",
                "message": {"source": source, "message": msg, **kwargs},
            }
        )


logger.add(handle_callback_log, filter=lambda record: record["level"].name == "INFO")
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "Info")
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "Warning")
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "Debug")
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "Error")
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "Exception")


# The Callback logger can hold a websocket that will allow the log message to be
# copied to the websocket as well as the logger
class CallbackLogger:
    def __init__(self, websocket: WebSocket):
        self.websocket = websocket
        self.logger = logger.bind(websocket=websocket)

    async def info(self, message, **kwargs):
        if kwargs:
            logger.bind(websocket=self.websocket, **kwargs).info(message)
        else:
            self.logger.info(message)

    async def warning(self, message, **kwargs):
        if kwargs:
            logger.bind(websocket=self.websocket, **kwargs).warning(message)
        else:
            self.logger.warning(message)

    async def debug(self, message, **kwargs):
        if kwargs:
            logger.bind(websocket=self.websocket, **kwargs).debug(message)
        else:
            self.logger.debug(message)

    async def error(self, message, **kwargs):
        if kwargs:
            logger.bind(websocket=self.websocket, **kwargs).error(message)
        else:
            self.logger.error(message)

    async def exception(self, message, **kwargs):
        if kwargs:
            logger.bind(websocket=self.websocket, **kwargs).error(message)
        else:
            self.logger.exception(message)

    def unbind(self):
        self.websocket = None
        self.logger = logger.bind()
