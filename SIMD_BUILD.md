# SIMD build selection

`configure` picks the SIMD mode at configure time. There is no runtime dispatch.

| `--enable-arch` | Target flags | Lanes (double) |
|---|---|---|
| `auto` (default) | probe picks `native` or `sse2` | 4 or 2 |
| `native` | `-march=native -ffp-contract=off` | 4 (AVX2 host) |
| `avx` | `-march=native -ffp-contract=off` | 4 |
| `sse2` | `-march=native -mno-avx -ffp-contract=off` | 2 |

- `auto` compiles a probe with `-march=native` and checks for `__AVX__`.
  AVX hosts get `native`, others get `sse2`. Float is never chosen
  automatically.
- `--enable-float` adds `-DCOSMO_FLOAT`. Opt-in only.
- The deprecated `--enable-sse2` and `--enable-avx` flags are read only
  when `--enable-arch` is not given. They map to the same modes.
- `-ffp-contract=off` is always set. It keeps `a*b+c` unfused, so output
  stays bit-exact with the fork, which uses no FMA.
- The probe prints `auto: host SIMD mode is ...` during configure. It warns
  on AVX-512 hosts, where the 8-lane path is not tested.

Examples:

```sh
./configure                       # auto
./configure --enable-arch=sse2    # override: SSE2 double
./configure --enable-arch=avx --enable-float
```

The wrapper `ChaNGa-test/Makefile` runs `auto` by default. Explicit
variants are listed in its header.
