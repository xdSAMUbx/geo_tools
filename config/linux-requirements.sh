#!/bin/bash
sudo apt update && sudo apt full-upgrade -y

sudo apt install mingw-w64 build-essential pkg-config \
    git curl -y