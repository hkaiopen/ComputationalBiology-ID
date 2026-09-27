#!/bin/bash
# Run all experiments sequentially.
# Results saved to results/ with timestamps.

set -e

echo "=== verify_grok_alignment ==="
python verify_grok_alignment.py

echo "=== exp_A ==="
python exp_A.py

echo "=== exp_B ==="
python exp_B.py

echo "=== exp_Bprime ==="
python exp_Bprime.py

echo "=== exp_C ==="
python exp_C.py

echo "=== exp_D ==="
python exp_D.py

echo "=== exp_Dprime ==="
python exp_Dprime.py

echo "=== exp_E ==="
python exp_E.py

echo "=== exp_F ==="
python exp_F.py

echo ""
echo "All experiments complete. Results in results/"
