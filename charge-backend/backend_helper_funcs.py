from loguru import logger
from fastapi import WebSocket
import asyncio
import json

RETROSYNTH_UNCONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise."
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
)

RETROSYNTH_CONSTRAINED_USER_PROMPT_TEMPLATE = (
    "Provide a retrosynthetic pathway for the target molecule {target_molecule}. "
    + "The pathway should be provided as a tuple of reactants as SMILES and the product as SMILES. "
    + "Perform only single step retrosynthesis. Make sure the SMILES strings are valid. "
    + "Use tools to verify the SMILES strings and diagnose any issues that arise. "
    + "The following reactant cannot be used in the retrosynthetic step: {constrained_reactant}. "
    + "Do the evaluation step-by-step. Propose a retrosynthetic step, then evaluate it. "
    + "If the evaluation fails, propose a new retrosynthetic step and evaluate it again. "
)


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
                            str_to_dict = json.loads(item.arguments)
                            if "smiles" in str_to_dict:
                                msg["smiles"] = str_to_dict["smiles"]
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
