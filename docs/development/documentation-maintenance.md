# Documentation maintenance

## Current authoritative files

InterMeta currently uses a mixed historical documentation layout. For this
release line, the tracked `NAMESPACE` file and tracked `man/*.Rd` files are the
authoritative package interface documentation and are maintained manually.
Roxygen comments present in some `R/*.R` files are useful source annotations,
but they are not yet the sole source of truth for the whole package.

Do not run a package-wide `roxygen2::roxygenise()` as a routine release step
until every exported function, S3 method, dataset, import, and manual detail has
been migrated and the resulting diff has been reviewed. A partial regeneration
could silently remove exports or overwrite richer hand-maintained help pages.

## Adding or changing an interface

For the current manual workflow:

1. update the function and its argument validation in `R/`;
2. add or change the matching `man/*.Rd` usage, arguments, return value,
   examples, and aliases;
3. update `NAMESPACE` for an export, S3 method, or import change;
4. update a runnable script in `inst/examples/` and the relevant short guide
   when user-facing behavior changes;
5. update `NEWS.md` for a release-visible change;
6. run documentation, installation, examples, tests, and package checks from
   the built source tarball.

For a dataset change, update its `.Rd` source and format sections, the public
generation script in `inst/scripts/`, its development mirror in `data-raw/`,
the relevant provenance object, fingerprints, and dataset tests together. The
two script copies must remain line-for-line content-equivalent; platform line
endings may differ. The two `update_prad_provenance.R` files are equivalent
thin wrappers; their shared, fail-closed implementation is maintained in
`inst/scripts/update_prad_provenance_impl.R` so it is also present in the
source tarball.

## Verification commands

From an R session at the package root, use:

```r
rd_files <- list.files("man", pattern = "[.]Rd$", full.names = TRUE)
rd_checks <- lapply(rd_files, tools::checkRd)
names(rd_checks) <- basename(rd_files)
Filter(length, rd_checks)
```

`checkRd()` reports diagnostics by severity. Review every result; existing
style diagnostics must not obscure a new syntax or cross-reference error.

Then build the source tarball and run `R CMD check` on that tarball. Inspect the
tarball contents to confirm that `inst/scripts/` reproduction scripts and
`inst/doc/` guides are present, and install it into a clean library to verify
that they appear under `system.file("scripts", package = "InterMeta")` and
`system.file("doc", package = "InterMeta")`.

## Future roxygen migration

If the project adopts roxygen as the single source of truth, make that a
dedicated migration: cover the entire exported API and all datasets, regenerate
`NAMESPACE` and `man/` with one pinned roxygen2 version, inspect the complete
diff, rerun all checks, and only then update `RoxygenNote`. Do not alternate
between manual and generated ownership after that migration.
