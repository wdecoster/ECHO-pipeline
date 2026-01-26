#!/bin/bash
set -euo pipefail

#----settings----
DB_NAME="echoDB_v1"
DEST_ROOT="resources"
DEST_DIR="${DEST_ROOT}/${DB_NAME}"
STAGING_DIR="${DEST_ROOT}/.${DB_NAME}.staging"
ZIP_PATH="${DEST_ROOT}/${DB_NAME}.zip"

# Download repeat catalog zip file from Zenodo to current directory
ZENODO="https://sandbox.zenodo.org/record/430049/files/${DB_NAME}.zip?download=1"


#----prep download----
mkdir -p "${DEST_ROOT}"

echo "[ECHO] Downloading ${DB_NAME}..."
wget -O "${ZIP_PATH}" "${ZENODO}"

echo "[ECHO] Extracting..."
rm -rf "${STAGING_DIR}"
mkdir -p "${STAGING_DIR}"
unzip -q -o "${ZIP_PATH}" -d "${STAGING_DIR}"

# We expect the zip to contain a top-level folder named echoDB_v1/
if [ ! -d "${STAGING_DIR}/${DB_NAME}" ]; then
  echo "[ECHO][ERROR] Expected top-level folder '${DB_NAME}/' inside the zip, but it was not found."
  echo "[ECHO][HINT] Ensure your zip contains '${DB_NAME}/TEs' and '${DB_NAME}/TRs' at its root."
  echo "[ECHO][HINT] Found these top-level entries in staging:"
  ls -lah "${STAGING_DIR}" || true
  exit 1
fi

echo "[ECHO] Installing to ${DEST_DIR}..."
rm -rf "${DEST_DIR}"
mv "${STAGING_DIR}/${DB_NAME}" "${DEST_DIR}"

rm -rf "${STAGING_DIR}"
rm -f "${ZIP_PATH}"

echo "[ECHO] Done. Installed catalogs in: ${DEST_DIR}"
echo "[ECHO] Quick check:"
echo "  - ${DEST_DIR}/TEs"
echo "  - ${DEST_DIR}/TRs"
