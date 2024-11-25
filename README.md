# An attack on Bloom-filter based (M)PSI schemes

This repository contains the proof-of-concept attack as described in "On the Insecurity of Bloom Filter-Based Private Set Intersections" by Jelle Vos, Jorrit van Assen, Tjitske Koster, Evangelia Anna Markatou, and Zekeriya Erkin, which can be found on the Cryptology ePrint archive as [2024/1901](https://ia.cr/2024/1901).
The attack demonstrates that many Bloom filter-based MPSI schemes underestimated the security-parameters.
The attacker exploits the knowledge about the hash functions, before the submission of its input set, to find elements that cause false positives with high probability.
This allows the attacker to learn elements from private input sets, that are not in the intersection.
For a detailed explanation of the attack and mitigations, please refer to the paper. 

## Overview

- `scr/` contains the Rust source code implementing the attack,
- `experiments/` contains code for executing the experiments. 

## To build and run

Run `cargo build --release` to compile the source code. Run `./target/release/bloom-goes-boom --help` to see the required parameters. Optionally run `python3 -u experiments/run_experiments.py` to replicate our results.
