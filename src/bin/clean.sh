#!/usr/bin/bash

echo "Try to clean GPU memory..."
nvidia-smi | grep python | awk '{print $5}' | xargs kill -9 || true

echo "Try to stop processes..."
ps -a | grep python | awk '{print $1}' | xargs kill -9 || true
ps -a | grep parallel | awk '{print $1}' | xargs kill -9 || true
ps -a | grep colabfold | awk '{print $1}' | xargs kill -9 || true

