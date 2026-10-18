# Syntax, dependency and package checks; does not load SEER data.
files <- c('run_analysis.R', 'check_code.R', list.files('work', pattern='\\.R$', full.names=TRUE))
for (f in files) parse(file=f, keep.source=FALSE)
required <- c('data.table', 'survival', 'splines', 'ggplot2', 'patchwork', 'svglite')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Required R packages unavailable: ', paste(missing, collapse=', '))
for (f in files) {
  code <- paste(readLines(f, warn=FALSE), collapse='\n')
  hits <- regmatches(code, gregexpr("work/[A-Za-z0-9_]+\\.R", code))[[1]]
  if (length(hits) && !all(file.exists(hits))) stop('Missing script dependency in ', f)
}
cat('PASS: ', length(files), ' R scripts parsed; packages and script dependencies available.\n', sep='')
