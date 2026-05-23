import base64
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from attachments import image_refs, validate_image_attachments


def image_attachment(**overrides):
    payload = {
        "id": "image-1",
        "name": "sample.png",
        "mimeType": "image/png",
        "sizeBytes": 4,
        "dataUrl": "data:image/png;base64," + base64.b64encode(b"fake").decode("ascii"),
        "createdAt": "2026-05-22T00:00:00Z",
    }
    payload.update(overrides)
    return payload


def test_validate_image_attachments_accepts_image_data_urls():
    attachments = validate_image_attachments({"attachments": [image_attachment()]})

    assert attachments[0]["id"] == "image-1"
    assert attachments[0]["mimeType"] == "image/png"
    assert attachments[0]["sizeBytes"] == 4


def test_validate_image_attachments_rejects_non_images():
    with pytest.raises(ValueError, match="must be an image"):
        validate_image_attachments(
            {
                "attachments": [
                    image_attachment(
                        mimeType="application/pdf",
                        dataUrl="data:application/pdf;base64,ZmFrZQ==",
                    )
                ]
            }
        )


def test_validate_image_attachments_rejects_more_than_five_files():
    with pytest.raises(ValueError, match="At most 5 images"):
        validate_image_attachments(
            {"attachments": [image_attachment(id=f"image-{i}") for i in range(6)]}
        )


def test_image_refs_omit_data_url():
    refs = image_refs([image_attachment()])

    assert refs == {
        "image-1": {
            "id": "image-1",
            "name": "sample.png",
            "mimeType": "image/png",
            "sizeBytes": 4,
        }
    }
