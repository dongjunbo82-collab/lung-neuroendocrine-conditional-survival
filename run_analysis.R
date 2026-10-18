# Run from the repository root; every stage gets its own R process.
core <- c('pnen_formal_OS.R', 'pnen_formal_competing.R', 'pnen_formal_sensitivity.R', 'pnen_formal_audit.R',
          'pnen_value_audit.R', 'pnen_clinical_standardization.R', 'pnen_clinical_competing.R')
all_steps <- c('pnen_formal_OS.R', 'pnen_horizon1.R', 'pnen_horizon_event_check.R',
 'pnen_formal_competing.R', 'pnen_formal_sensitivity.R', 'pnen_formal_descriptive.R',
 'pnen_formal_audit.R', 'pnen_sensitivity_contrasts.R', 'pnen_reporting_descriptives.R',
 'pnen_figures.R', 'pnen_extension_audit.R', 'pnen_age_support_audit.R',
 'pnen_extension_descriptive.R', 'pnen_extension_models.R', 'pnen_extension_figures.R',
 'pnen_lung_additions.R', 'pnen_lung_figures.R', 'pnen_value_audit.R',
 'pnen_clinical_standardization.R', 'pnen_clinical_competing.R', 'pnen_tlcr_figure.R')
args <- commandArgs(trailingOnly=TRUE)
if (identical(args, '--list')) {
  cat('Core:\n', paste(core, collapse='\n'), '\n\nAll:\n', paste(all_steps, collapse='\n'), '\n', sep='')
} else {
  if (length(args)!=1L || !args %in% c('--core', '--all')) stop('Use --list, --core or --all')
  if (!file.exists('work/pnen_formal_OS.R')) stop('Run from the repository root')
  project <- normalizePath(getwd(), winslash='/')
  if (nzchar(Sys.getenv('PNEN_PROJECT')) && normalizePath(Sys.getenv('PNEN_PROJECT'), winslash='/') != project)
    stop('PNEN_PROJECT must refer to this repository root')
  for (key in c('SEER_INPUT', 'SEER_EXTENDED'))
    if (!nzchar(Sys.getenv(key)) || !file.exists(Sys.getenv(key))) stop('Set ', key, ' to an accessible SEER export')
  source('check_code.R', local=TRUE)
  Sys.setenv(PNEN_PROJECT=project)
  dir.create('outputs/pnen_manuscript', recursive=TRUE, showWarnings=FALSE)
  writeLines(capture.output(sessionInfo()), 'outputs/session_info.txt')
  steps <- if (args=='--core') core else all_steps
  exe <- file.path(R.home('bin'), if (.Platform$OS.type=='windows') 'Rscript.exe' else 'Rscript')
  for (step in steps) {
    if (step=='pnen_lung_figures.R') {
      # The retained plotting script expects this aggregate-file alias.
      copied <- file.copy('outputs/pnen_extension/extension_estimates.csv',
                          'outputs/pnen_manuscript/Supplement_extension_extension_estimates.csv', overwrite=TRUE)
      if (!copied) stop('Cannot prepare the extension-estimate input for plotting')
    }
    cat('Running ', step, '\n', sep='')
    status <- system2(exe, args=shQuote(file.path('work', step)))
    if (status!=0L) stop('Analysis failed: ', step)
  }
}
