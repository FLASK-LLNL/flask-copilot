import sys
from fastapi import WebSocket
from loguru import logger

# Define the callback function - will send message to the websocket if it is provided
async def handle_callback_log(message):
    record = message.record
    websocket = record["extra"].get("websocket", None)
    if websocket:
        await websocket.send_json(
            {
                "type": "response",
                "message": message,
            }
        )
    sys.stdout.write(message)
        
logger.add(handle_callback_log, filter=lambda record: record["level"].name == "INFO")

# The Callback logger can hold a websocket that will allow the log message to be
# copied to the websocket as well as the logger
class callback_logger:
     def __init__(self, websocket: WebSocket):
         self.websocket = websocket
         self.logger = logger.bind(websocket=websocket)
     def info(self, message):
         self.logger.info(message)

