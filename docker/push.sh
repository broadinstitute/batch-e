#!/bin/bash
set -ex

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

SHA="$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
if ! git -C "$REPO_ROOT" diff --quiet HEAD -- .; then
    SHA="${SHA}-dirty"
fi

REPO="us-docker.pkg.dev/broad-dsde-methods/batch-e/batch-e"

docker buildx build \
    -t "${REPO}:${SHA}" \
    -t "${REPO}:latest" \
    --platform linux/amd64 \
    --build-arg "GIT_SHA=${SHA}" \
    --push \
    -f "$SCRIPT_DIR/Dockerfile" \
    "$REPO_ROOT"
