"""Resumably download and validate Cell Model Passports RNA-seq ZIP files."""

from __future__ import annotations

import csv
import time
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw" / "rnaseq"
PROCESSED = ROOT / "processed"
STATUS_PATH = PROCESSED / "download_status.csv"


def valid_zip(path: Path, model_id: str) -> bool:
    if not path.exists() or path.stat().st_size < 1000:
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:
                return False
            names = archive.namelist()
            if names != [f"{model_id}_rnaseq.csv"]:
                return False
            with archive.open(names[0]) as stream:
                header = stream.readline().decode("utf-8").strip().split(",")
            return {"model_id", "gene_symbol", "rsem_tpm"}.issubset(header)
    except (OSError, zipfile.BadZipFile, UnicodeDecodeError):
        return False


def download_one(record: dict[str, str]) -> dict[str, str | int]:
    model_id = record["model_id"]
    url = record["rnaseq_url"]
    destination = RAW / f"{model_id}_rnaseq.csv.zip"
    partial = destination.with_suffix(destination.suffix + ".part")
    if valid_zip(destination, model_id):
        return {"model_id": model_id, "status": "existing", "bytes": destination.stat().st_size, "message": ""}

    last_error = ""
    for attempt in range(1, 5):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "STAT293-GDSC-pilot/1.0"})
            with urllib.request.urlopen(request, timeout=90) as response, partial.open("wb") as out:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    out.write(block)
            partial.replace(destination)
            if not valid_zip(destination, model_id):
                raise RuntimeError("downloaded file failed ZIP/content validation")
            return {"model_id": model_id, "status": "downloaded", "bytes": destination.stat().st_size, "message": ""}
        except (OSError, RuntimeError, urllib.error.URLError) as error:
            last_error = f"attempt {attempt}: {error}"
            if partial.exists():
                partial.unlink()
            time.sleep(min(2**attempt, 10))
    return {"model_id": model_id, "status": "failed", "bytes": 0, "message": last_error}


def write_status(rows: list[dict[str, str | int]]) -> None:
    rows = sorted(rows, key=lambda row: str(row["model_id"]))
    with STATUS_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["model_id", "status", "bytes", "message"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(PROCESSED / "model_manifest.csv", dtype=str)
    records = manifest[["model_id", "rnaseq_url"]].to_dict("records")
    completed: list[dict[str, str | int]] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(download_one, record): record for record in records}
        for number, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            completed.append(result)
            write_status(completed)
            print(f"[{number:03d}/{len(records):03d}] {result['model_id']}: {result['status']}", flush=True)
    failures = [row for row in completed if row["status"] == "failed"]
    if failures:
        raise RuntimeError(f"{len(failures)} RNA-seq downloads failed; see {STATUS_PATH}")


if __name__ == "__main__":
    main()
