# Fork-specific helper: list genome assembly builds with bundled cytoband data

## Output:
### character vector of assembly names with a matching data.frame in R/sysdata.rda

## Required by:
### pcf
### multipcf
### aspcf
### winsorize

## Requires:
### none

validAssemblies <- function() {
  ns <- parent.env(environment())
  objs <- ls(ns)
  objs[vapply(objs, function(o) is.data.frame(get(o, envir = ns)), logical(1))]
}
