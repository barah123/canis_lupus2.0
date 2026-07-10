# Contributing to CanisLupus 2.0

Thanks for your interest in contributing to CanisLupus 2.0! This document
covers how to report bugs, request features, and submit changes.

## Reporting bugs

If you find a bug, please open an issue on GitHub and include:

- A short description of what you expected to happen vs. what actually happened
- Steps to reproduce (ideally with a small example dataset, or one of the
  sample datasets bundled with the app)
- Your R version (`R.version.string`) and operating system
- Any error messages or console output

## Requesting features

Feature requests are welcome via GitHub issues. Please describe the use case
(what analysis or workflow you're trying to accomplish) so we can evaluate
fit with the project's scope.

## Submitting changes

1. Fork the repository and create a new branch for your change
   (`git checkout -b fix/short-description`).
2. Make your changes. Please keep pull requests focused on a single fix or
   feature where possible.
3. If your change affects app behavior, please run the existing tests
   locally before submitting:
   ```r
   testthat::test_dir("tests/testthat")
   ```
   and add a test if you're fixing a bug or adding functionality.
4. Update the README if your change affects installation, usage, or
   supported data formats.
5. Open a pull request describing what the change does and why.

## Code style

- Follow the [tidyverse style guide](https://style.tidyverse.org/) for R code
  where practical.
- Keep functions focused and commented, especially around data
  transformations (e.g., ASV table filtering, distance matrix computation)
  where assumptions about input format matter.

## Getting help

For questions about using the app (as opposed to bug reports), please open a
GitHub issue with the `question` label rather than emailing directly, so
answers are visible to other users with the same question.

## Code of conduct

Participants are expected to treat each other with respect. Harassment or
discriminatory language/behavior of any kind will not be tolerated. Instances
of unacceptable behavior can be reported by opening an issue or contacting
the maintainer directly.
