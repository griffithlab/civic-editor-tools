#!/usr/bin/env bash
set -euo pipefail

# Self-locate the repo root so this script works regardless of the caller's cwd, and regardless
# of install/docker/ moving again in the future.
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

IMAGE="griffithlab/civic-editor-tools"
VERSION=$(cat "$REPO_ROOT/RELEASE" | tr -d '[:space:]')

# Verify Docker Hub authentication before wasting time on the build
CREDS_STORE=$(python3 -c "import json; d=json.load(open('$HOME/.docker/config.json')); print(d.get('credsStore',''))" 2>/dev/null || true)
if [[ -n "$CREDS_STORE" ]]; then
    AUTH_CHECK=$(echo "https://index.docker.io/v1/" | docker-credential-"${CREDS_STORE}" get 2>/dev/null || true)
else
    AUTH_CHECK=$(python3 -c "import json; d=json.load(open('$HOME/.docker/config.json')); print(d.get('auths',{}).get('https://index.docker.io/v1/',{}).get('auth',''))" 2>/dev/null || true)
fi
if [[ -z "$AUTH_CHECK" ]]; then
    echo "Error: not logged in to Docker Hub. Run 'docker login' first." >&2
    exit 1
fi

echo "Building and pushing ${IMAGE}:${VERSION} and ${IMAGE}:latest"
echo "Platforms: linux/amd64, linux/arm64"

# Ensure a buildx builder with multi-platform support is active
if ! docker buildx inspect multiplatform-builder &>/dev/null; then
    docker buildx create --name multiplatform-builder --use
else
    docker buildx use multiplatform-builder
fi

docker buildx build \
    --platform linux/amd64,linux/arm64 \
    --tag "${IMAGE}:${VERSION}" \
    --tag "${IMAGE}:latest" \
    --push \
    -f "$SCRIPT_DIR/Dockerfile" \
    "$REPO_ROOT"

echo "Done: ${IMAGE}:${VERSION} and ${IMAGE}:latest pushed to Docker Hub"
