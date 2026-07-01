# Filematrix backend design note

This note records the current filematrix backend strategy for bigPLSR and
bigPCAcpp. The goal is benchmark validation of row-block file-backed workflows,
not replacement of the existing bigmemory backends.

## 1. Problem statement

bigmemory remains supported and valuable. It provides compatibility with the
existing C++ paths and is appropriate for many local, file-backed, and
random-access workflows.

The stress campaign showed a separate storage-layer issue: mmap-backed writes
to very large file-backed matrices can be fragile or slow on shared
filesystems. Repeated failures were observed during benchmark data
materialization, with SIGBUS / signal 7 in SetMatrixRows or SetMatrixElements.
Those failures occurred while constructing or assigning file-backed benchmark
inputs, not inside the corrected streaming PCA or PLS algorithm logic.

These events should be classified as storage/infrastructure failures unless the
same signal occurs inside the PCA or PLS algorithm code itself. When input
materialization was controlled, the corrected streaming algorithms continued to
work, including streaming NIPALS PLS1/PLS2 at p = 100000 predictors.

## 2. Storage-backend distinction

```text
bigmemory:
  memory-mapped file-backed matrices;
  existing C++ algorithm compatibility;
  useful for existing workflows and some random-access workloads;
  can depend on OS page-cache and mmap semantics.

filematrix:
  explicit file-backed matrix access through R-level reads/writes;
  better aligned with row-block streaming experiments;
  avoids using bigmemory mmap for filematrix input reads;
  not claimed to be universally faster.
```

Benchmark reports should keep three labels separate:

```text
generator:
  the mechanism used to create the input files.

storage_backend:
  the on-disk matrix representation used for generated data.

algorithm_input_backend:
  the backend consumed by the fitting or projection algorithm.
```

Separating these labels prevents a generation or attachment failure from being
misreported as an algorithm failure.

## 3. Implemented backend strategy

### bigPLSR

- Existing bigmemory backends remain in place.
- An optional filematrix row-block provider is available.
- filematrix remains in Suggests, not Imports.
- pls_fit() supports backend = "filematrix" experimentally for NIPALS PLS1 and
  PLS2 only.
- The filematrix NIPALS implementation reads row blocks through the provider and
  avoids bigmemory mmap input access.
- SIMPLS, kernel PLS, wide kernel PLS, RKHS, and KF-PLS remain outside the
  filematrix backend scope.
- Existing bigmemory code paths are not changed by the filematrix backend.

### bigPCAcpp

- Existing bigmemory backends remain in place.
- An optional filematrix row-block provider is available.
- filematrix remains in Suggests, not Imports.
- pca_spca_stream_filematrix() implements streaming SPCA for filematrix inputs.
- pca_scores_stream_filematrix() implements dense score projection by row-block
  reads.
- pca_stream_filematrix() implements exact covariance PCA for moderate-p
  filematrix benchmark validation.
- pca_stream_filematrix() forms a p x p covariance matrix and is guarded by
  getOption("bigPCAcpp.filematrix.max_cov_gb", 2).
- Very-wide filematrix PCA should use pca_spca_stream_filematrix(), not exact
  covariance PCA.
- Existing bigmemory code paths are not changed by the filematrix backend.

## 4. Native C++ integration decision

Native C++ integration is deferred.

The first implementation intentionally uses R-level filematrix row-block reads
and R-level orchestration. Native C++ support should be considered only if
R-level filematrix benchmarks show material reliability or performance benefits
and if the file format/API can be accessed safely from C++ without compromising
CRAN portability.

Native C++ integration is not required to validate the main storage-layer
hypothesis. The first question is whether replacing mmap-backed input
materialization and reads by explicit row-block filematrix reads improves
reliability on shared filesystems.

## 5. Benchmark plan

For the benchmark layer:

- Compare bigmemory and filematrix generators on tiny_smoke, small_safe, and
  wide_quota_test.
- Keep explicit generator, storage_backend, and algorithm_input_backend labels.

For bigPLSR:

- Compare filematrix NIPALS PLS1/PLS2 against dense or bigmemory references on
  small data.
- Track coefficients, prediction agreement, elapsed time, status, scratch
  footprint, and failure mode.

For bigPCAcpp:

- Compare pca_spca_stream_filematrix() against dense pca_spca_R() on small data.
- Compare pca_scores_stream_filematrix() against dense score projection.
- Compare pca_stream_filematrix() against prcomp() or existing exact PCA on
  moderate-p data.
- Do not run exact filematrix PCA at very-wide p unless the covariance-size
  guard is deliberately raised and sufficient memory is available.

General metrics:

- elapsed time;
- user/system time;
- scratch footprint;
- status;
- failure mode;
- numerical parity against reference;
- whether failure occurred in generation, storage attachment, algorithm fitting,
  or score projection.

## 6. Failure classification

A SIGBUS / signal 7 occurring during file-backed matrix generation or assignment
is classified as a storage/materialization failure.

A SIGBUS / signal 7 occurring inside a streaming algorithm kernel is classified
separately as an algorithm/runtime failure.

A successful filematrix generator run followed by successful algorithm execution
supports the diagnosis that earlier failures came from mmap-backed
materialization rather than from the PCA/PLS streaming algorithms.

## 7. Scope limitations

filematrix support is experimental in these packages. It is not claimed to be
universally faster than bigmemory.

The exact PCA filematrix path is for moderate-p benchmark validation only. The
SPCA and PLS NIPALS filematrix paths are the relevant very-wide streaming paths.

Native C++ integration is intentionally deferred.
