from typing import Any
import base64
import binascii


MAX_IMAGE_ATTACHMENTS = 5
MAX_IMAGE_ATTACHMENT_BYTES = 5 * 1024 * 1024


def validate_image_attachments(data: dict[str, Any]) -> list[dict[str, Any]]:
    attachments = data.get("attachments") or []
    if not attachments:
        return []
    if not isinstance(attachments, list):
        raise ValueError("attachments must be a list")
    if len(attachments) > MAX_IMAGE_ATTACHMENTS:
        raise ValueError(f"At most {MAX_IMAGE_ATTACHMENTS} images can be attached")

    normalized: list[dict[str, Any]] = []
    for index, attachment in enumerate(attachments):
        if not isinstance(attachment, dict):
            raise ValueError(f"Attachment {index + 1} must be an object")

        data_url = attachment.get("dataUrl")
        mime_type = attachment.get("mimeType")
        if not isinstance(data_url, str) or not data_url.startswith("data:"):
            raise ValueError(f"Attachment {index + 1} must be a data URL")
        if not isinstance(mime_type, str) or not mime_type.startswith("image/"):
            raise ValueError(f"Attachment {index + 1} must be an image")

        header, separator, encoded = data_url.partition(",")
        if separator != "," or ";base64" not in header:
            raise ValueError(f"Attachment {index + 1} must be base64 encoded")

        header_mime_type = header.removeprefix("data:").split(";", 1)[0]
        if header_mime_type != mime_type:
            raise ValueError(f"Attachment {index + 1} has inconsistent MIME metadata")

        try:
            decoded_size = len(base64.b64decode(encoded, validate=True))
        except (binascii.Error, ValueError) as exc:
            raise ValueError(f"Attachment {index + 1} is not valid base64") from exc
        if decoded_size > MAX_IMAGE_ATTACHMENT_BYTES:
            raise ValueError(f"Attachment {index + 1} exceeds 5 MB")

        normalized.append(
            {
                "id": str(attachment.get("id") or f"image_{index + 1}"),
                "name": str(attachment.get("name") or f"Image {index + 1}"),
                "mimeType": mime_type,
                "sizeBytes": decoded_size,
                "dataUrl": data_url,
                "createdAt": str(attachment.get("createdAt") or ""),
            }
        )

    return normalized


def image_refs(attachments: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        attachment["id"]: {
            "id": attachment["id"],
            "name": attachment["name"],
            "mimeType": attachment["mimeType"],
            "sizeBytes": attachment["sizeBytes"],
        }
        for attachment in attachments
    }
