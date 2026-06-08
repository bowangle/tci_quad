#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EIGEN_DIR="$SCRIPT_DIR/eigen"
JSON_DIR="$SCRIPT_DIR/nlohmann_json"
SPDLOG_DIR="$SCRIPT_DIR/spdlog"

echo "==> Cloning Eigen..."
git clone https://gitlab.com/libeigen/eigen.git "$EIGEN_DIR"

echo "==> Downloading nlohmann/json..."
mkdir -p "$JSON_DIR/nlohmann"

curl -L https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp \
    -o "$JSON_DIR/nlohmann/json.hpp"

echo "==> Cloning spdlog..."
if [ ! -d "$SPDLOG_DIR" ]; then
    git clone https://github.com/gabime/spdlog.git "$SPDLOG_DIR"
else
    echo "spdlog already exists, skipping."
fi

echo "==> Done!"

echo ""
echo "Eigen cloned to:"
echo "$EIGEN_DIR"

echo ""
echo "nlohmann/json installed to:"
echo "$JSON_DIR"

echo ""
echo "Use it by compiling with:"
echo "  -I $EIGEN_DIR"
echo "  -I $JSON_DIR"
echo ""
echo "And in code:"
echo "  #include <nlohmann/json.hpp>"