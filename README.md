# QGDecoder

**QGDecoder** is a graph-based bounded-distance quantum decoder supporting arbitrary stabilizer codes
- An $[[N,k,d]]$ quantum code corrects upto weight $t=(d-1)/2$ errors with certainty.   
- QGDecoder takes in a tunable target weight $T$. If $T\leq t$ it ensures all errors with weight $w\leq T$ is corrected.
- Graph state based bounded distance decoding.
- Additive.h supports both CSS and non-CSS codes.
- Use CSS.h for improved efficiency with CSS codes.

## Package features and requirements
- Header only package. Make sure the parent folder QGDecoder is accessible. Enjoy decoding!
- Written in C++ (use C++17 preferably).
- Linear algebra powered by Armadillo. Install from [here](http://arma.sourceforge.net/) to use QGDecoder.

## Available quantum codes
- non-CSS [optimal codes](https://www.codetables.de/) of distance $d=3,5,7,9,11$.
- non-CSS [XZZX code](https://www.nature.com/articles/s41467-021-22274-1) on $d\times d$ square lattice.
- CSS triangular color code on triangular lattice of length $d$.
- CSS rotated surface codes on $d\times d$ square lattice.  

## Attribution

When using QGDecoder, please cite our [paper](https://arxiv.org/abs/):

```bibtex
@article{Harikrishnan2026QGDecoder,
  title         = {A graph-aware bounded-distance decoder for all stabilizer codes},
  author        = {Harikrishnan K J and Amit Kumar Pal},
  journal       = {arXiv preprint},
  year          = {2026},
  eprint        = {},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
  url           = {https://arxiv.org/abs}
}
