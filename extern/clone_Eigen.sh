#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EIGEN_DIR="$SCRIPT_DIR/eigen"

echo "==> Cloning Eigen..."
git clone https://gitlab.com/libeigen/eigen.git "$EIGEN_DIR"

echo "==> Done!"
echo "Eigen cloned to:"
echo "$EIGEN_DIR"

echo ""
echo "Use it by compiling with:"
echo "  -I $EIGEN_DIR"
echo ""
echo "And in code:"
echo "  #include <Eigen/Dense>"
